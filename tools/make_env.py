#!/usr/bin/env python3
"""Legacy reproduction only: use a separate baseline checkout, never overwrite frozen fixtures.

Generate a task ``environment/`` from ``solver-python`` plus an ablation spec.

**``environment/`` is generated, never hand-edited.**  Each task declares an
``ablation.yaml`` saying what to remove and how to rewrite what is left; this
applies it to the current solver.  When the solver changes, every environment
is regenerated with one command, and a patch that no longer matches fails
loudly instead of drifting silently.

The scrub is an *assertion*, not a text substitution.  Mangling the word
"acceleration" out of a comment leaves a comment that reads oddly and still
tells the agent where to look; the build instead refuses to finish while any
forbidden term survives, and names every file and line, so the ablation gets
fixed properly.

Usage::

    python tools/make_env.py tasks/t01-square [--check]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

import yaml

SKIP_DIRS = {"__pycache__", ".git", ".pytest_cache", "out", "pybte.egg-info",
             ".ruff_cache", ".mypy_cache"}


def _iter_files(root: Path):
    for p in sorted(root.rglob("*")):
        if p.is_file() and not any(d in p.parts for d in SKIP_DIRS):
            yield p


def copy_source(src: Path, dst: Path) -> None:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns(*SKIP_DIRS))


def remove_paths(env: Path, paths) -> None:
    for rel in paths or []:
        for target in sorted(env.glob(rel)):
            if target.is_dir():
                shutil.rmtree(target)
            elif target.exists():
                target.unlink()


def apply_patches(env: Path, patches) -> None:
    for i, patch in enumerate(patches or []):
        f = env / patch["file"]
        if not f.exists():
            raise SystemExit(f"patch {i}: {patch['file']} does not exist")
        text = f.read_text()
        find, repl = patch["find"], patch.get("replace", "")
        n = text.count(find)
        want = patch.get("count", 1)
        if n != want:
            raise SystemExit(
                f"patch {i} on {patch['file']}: expected {want} match(es), found {n}.\n"
                f"The solver has moved under the ablation -- fix ablation.yaml.\n"
                f"--- find ---\n{find[:400]}")
        f.write_text(text.replace(find, repl))


def truncate(env: Path, specs) -> None:
    """Cut a file at the first line of ``marker`` -- for lopping a trailing
    block (a whole group of functions) without pasting it into the spec."""
    for spec in specs or []:
        f = env / spec["file"]
        text = f.read_text()
        i = text.find(spec["marker"])
        if i < 0:
            raise SystemExit(f"truncate: marker not found in {spec['file']}")
        f.write_text(text[:i].rstrip() + "\n")


def write_new_files(env: Path, files) -> None:
    for rel, body in (files or {}).items():
        p = env / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(body)


def scrub_check(env: Path, terms, exempt) -> list[str]:
    """Every forbidden term, in contents and in filenames.  Returns violations."""
    if not terms:
        return []
    pat = re.compile("|".join(rf"\b{re.escape(t)}\b" for t in terms), re.I)
    exempt = {str(e) for e in (exempt or [])}
    bad = []
    for p in _iter_files(env):
        rel = str(p.relative_to(env))
        if rel in exempt:
            continue
        if pat.search(p.name):
            bad.append(f"{rel}: forbidden term in FILENAME")
        try:
            lines = p.read_text().splitlines()
        except (UnicodeDecodeError, OSError):
            continue
        for n, line in enumerate(lines, 1):
            m = pat.search(line)
            if m:
                bad.append(f"{rel}:{n}: {m.group(0)!r}  |  {line.strip()[:90]}")
    return bad


def manifest(env: Path) -> dict:
    out = {}
    for p in _iter_files(env):
        out[str(p.relative_to(env))] = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", type=Path)
    ap.add_argument("--check", action="store_true",
                    help="only report scrub violations for an existing environment")
    args = ap.parse_args(argv)

    task = args.task.resolve()
    spec = yaml.safe_load((task / "ablation.yaml").read_text())
    env = task / "environment"

    if not args.check:
        src = (task / spec["source"]).resolve()
        if not src.is_dir():
            raise SystemExit(f"source {src} not found")
        print(f"source      {src}")
        copy_source(src, env)
        remove_paths(env, spec.get("remove_paths"))
        apply_patches(env, spec.get("patch"))
        truncate(env, spec.get("truncate"))
        write_new_files(env, spec.get("new_files"))
        if spec.get("canary"):
            (env / "CANARY").write_text(spec["canary"] + "\n")

    bad = scrub_check(env, spec.get("scrub_terms"), spec.get("scrub_exempt"))
    n_files = sum(1 for _ in _iter_files(env))
    if bad:
        print(f"\nSCRUB FAILED -- {len(bad)} violation(s):", file=sys.stderr)
        for b in bad[:60]:
            print("  " + b, file=sys.stderr)
        if len(bad) > 60:
            print(f"  ... and {len(bad) - 60} more", file=sys.stderr)
        return 1

    if not args.check:
        (task / "manifest.json").write_text(json.dumps(manifest(env), indent=2))
    print(f"environment {env}  ({n_files} files)   scrub clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
