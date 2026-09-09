#!/usr/bin/env python3
"""Build a task's oracle solution and emit it as a patch against ``environment/``.

The oracle has to clear the verifier five times out of five.  If it does not,
the task is unsolvable as specified and no agent result from it means
anything -- which is why this runs before any rollout, not after.

The solution tree is built from the same source the environment is, with the
ablation's *shaping* applied (case files, layout) but not its *removals*, so
the two differ exactly by the capability the task asks for.

Usage::

    python tools/make_oracle.py tasks/t01-diffusive-acceleration
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_env import (SKIP_DIRS, apply_patches, copy_source, remove_paths,  # noqa: E402
                      write_new_files)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", type=Path)
    args = ap.parse_args(argv)

    task = args.task.resolve()
    spec = yaml.safe_load((task / "ablation.yaml").read_text())
    osp = yaml.safe_load((task / "oracle" / "solution.yaml").read_text())
    src = (task / spec["source"]).resolve()
    env = task / "environment"
    sol = task / "oracle" / "solution"

    copy_source(src, sol)
    remove_paths(sol, osp.get("remove_paths"))
    apply_patches(sol, osp.get("patch"))
    write_new_files(sol, osp.get("new_files"))

    # the environment's tests carry the same shaping; keep them identical so
    # the patch shows only what the task is actually about
    if (env / "tests").is_dir():
        shutil.rmtree(sol / "tests", ignore_errors=True)
        shutil.copytree(env / "tests", sol / "tests",
                        ignore=shutil.ignore_patterns(*SKIP_DIRS))
    for extra in ("CANARY",):
        if (env / extra).exists():
            shutil.copy2(env / extra, sol / extra)

    diff = subprocess.run(
        ["diff", "-ruN", "-x", "__pycache__", "-x", "*.pyc",
         str(env.relative_to(task)), str(sol.relative_to(task))],
        cwd=task, capture_output=True, text=True)
    if diff.returncode not in (0, 1):
        print(diff.stderr, file=sys.stderr)
        return 2
    out = task / "oracle" / "patch.diff"
    out.write_text(diff.stdout)
    n_files = len([l for l in diff.stdout.splitlines() if l.startswith("diff -ruN")])
    print(f"solution  {sol}")
    print(f"patch     {out}  ({n_files} files, "
          f"{len(diff.stdout.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
