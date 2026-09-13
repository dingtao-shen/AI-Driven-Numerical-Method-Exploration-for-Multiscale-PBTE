#!/usr/bin/env python3
"""Re-grade archived rollouts under the verifier as it is now.

``run_rollouts.py`` archives every solution as a tarball of the environment
the agent left behind.  When the verifier changes -- a rule refined, a
grading fault fixed -- the archived solutions are scored again here, so the
numbers of a run stay comparable without spending rollouts.  The agent is
never re-run; only the grading is.

Each tarball is unpacked into a scratch directory outside the repository
(the same rule as for rollouts), graded by ``verifier/checks.py`` exactly as
a fresh trial would be, and its ``trial_XX_score.json`` and row in
``summary.json`` replaced.  ``--oracle`` re-scores the task's oracle into
``oracle_score.json`` alongside.

Usage::

    python tools/rescore_rollouts.py experiments/results/t01-square-F3 \\
        --task tasks/t01-square --sandbox-root ~/rollout-sandboxes/rescore \\
        --trials 1 2 3 --oracle
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def grade(task: Path, env: Path, out_json: Path) -> dict:
    r = subprocess.run(
        [sys.executable, str(task / "verifier" / "checks.py"),
         "--env", str(env), "--json", str(out_json)],
        capture_output=True, text=True)
    if not out_json.exists():
        return {"gate": False, "score": None, "families": {},
                "verifier_error": (r.stderr or r.stdout)[-2000:]}
    return json.loads(out_json.read_text())


def row_from(res: dict) -> dict:
    cells = [c for f in res.get("families", {}).values() for c in f["cells"]]
    return {"gate": bool(res.get("gate")), "score": res.get("score"),
            "families": {k: (f["gate"], f["score"]) for k, f in res.get("families", {}).items()},
            "cells_passed": sum(1 for c in cells if c.get("passed")), "cells": len(cells)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", type=Path)
    ap.add_argument("--task", type=Path, required=True)
    ap.add_argument("--sandbox-root", type=Path, required=True,
                    help="scratch directory for the unpacked solutions; outside the repository")
    ap.add_argument("--trials", type=int, nargs="*", default=None,
                    help="trial numbers to re-score (default: every archived trial)")
    ap.add_argument("--oracle", action="store_true", help="re-score the oracle too")
    args = ap.parse_args(argv)

    task, results = args.task.resolve(), args.results.resolve()
    root = args.sandbox_root.expanduser().resolve()
    if ROOT in root.parents or root == ROOT:
        print("refusing to unpack solutions inside the repository", file=sys.stderr)
        return 2
    root.mkdir(parents=True, exist_ok=True)

    trials = args.trials or sorted(
        int(p.name[6:8]) for p in results.glob("trial_*_environment.tar.gz"))
    summary_path = results / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {"trials": []}
    rows = {r["trial"]: r for r in summary["trials"]}

    if args.oracle:
        res = grade(task, task / "oracle" / "solution", results / "oracle_score.json")
        r = row_from(res)
        print(f"oracle:   gate {'PASS' if r['gate'] else 'FAIL'}   "
              f"score {r['score']:.1f}x   {r['cells_passed']}/{r['cells']} cells"
              if r["score"] else f"oracle:   gate FAIL   {r['cells_passed']}/{r['cells']} cells",
              flush=True)

    for i in trials:
        tb = results / f"trial_{i:02d}_environment.tar.gz"
        if not tb.exists():
            print(f"trial {i}: no archive", flush=True)
            continue
        box = root / f"trial_{i:02d}"
        if box.exists():
            shutil.rmtree(box)
        box.mkdir(parents=True)
        with tarfile.open(tb) as tf:
            tf.extractall(box, filter="data")
        res = grade(task, box / "environment", results / f"trial_{i:02d}_score.json")
        r = row_from(res)
        prev = rows.get(i, {"trial": i})
        prev_score = prev.get("score")
        rows[i] = {**prev, **r, "rescored": True}
        sc = f"{r['score']:.1f}x" if r["score"] else "-"
        was = f"{prev_score:.1f}x" if prev_score else "-"
        print(f"trial {i}: gate {'PASS' if r['gate'] else 'FAIL'}   score {sc} (was {was})   "
              f"{r['cells_passed']}/{r['cells']} cells", flush=True)

    summary["trials"] = [rows[k] for k in sorted(rows)]
    counted = [r for r in summary["trials"] if not r.get("voided")]
    summary["solved"] = sum(1 for r in counted if r["gate"])
    summary["voided"] = len(summary["trials"]) - len(counted)
    summary_path.write_text(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
