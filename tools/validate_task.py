#!/usr/bin/env python3
"""Gate a task before any rollout is spent on it.

Five things, all of which have to hold:

1. **scrub clean** -- no forbidden term survives in the environment.
2. **oracle scores 1.0**, k times out of k.  If it does not, the task is
   unsolvable as specified and no agent result from it means anything.
3. **null baseline scores 0.0**, k times out of k.  If it does not, the task
   is already solved by the code it ships.
4. **the shrink preserved the pathology** -- the solver as shipped must still
   degrade by orders of magnitude across the stiffness ladder.  A shrink that
   accidentally lets the baseline converge quickly destroys the task, and it
   destroys it silently.
5. **verification fits the wall-clock budget**, with room for an agent to
   have run it many times while iterating.

k > 1 is not about statistics -- the solver is deterministic.  It is about
catching nondeterminism that should not be there: thread races, hash
ordering, anything that makes a score depend on the run.

Usage::

    python tools/validate_task.py tasks/t01-diffusive-acceleration [--k 5]
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def score(task: Path, env: Path, tag: str):
    with tempfile.TemporaryDirectory() as td:
        js = Path(td) / "r.json"
        t0 = time.perf_counter()
        subprocess.run([sys.executable, str(task / "verifier" / "checks.py"),
                        "--env", str(env), "--json", str(js)],
                       capture_output=True, text=True)
        wall = time.perf_counter() - t0
        if not js.exists():
            return None, wall, None
        d = json.loads(js.read_text())
        cells = [c for f in d["families"].values() for c in f["cells"]]
        return bool(d["gate"]), wall, cells


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", type=Path)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--budget-seconds", type=float, default=900.0)
    args = ap.parse_args(argv)
    task = args.task.resolve()
    ok = True

    print("1. scrub")
    r = subprocess.run([sys.executable, str(ROOT / "tools" / "make_env.py"),
                        str(task), "--check"], capture_output=True, text=True)
    print("   " + (r.stdout.strip().splitlines() or ["(no output)"])[-1])
    if r.returncode != 0:
        print(r.stderr[:2000], file=sys.stderr)
        ok = False

    walls = []
    for tag, env, want in (("2. oracle", task / "oracle" / "solution", True),
                           ("3. null baseline", task / "environment", False)):
        print(f"\n{tag}   (want gate {'PASS' if want else 'FAIL'} x{args.k})")
        got = []
        for i in range(args.k):
            s, wall, cases = score(task, env, tag)
            walls.append(wall)
            got.append(s)
            n = sum(1 for c in (cases or []) if c.get("passed"))
            print(f"   run {i+1}: gate {'PASS' if s else 'FAIL'}   "
                  f"{n}/{len(cases or [])} cells   {wall:.0f}s", flush=True)
        if any(s != want for s in got):
            print(f"   FAILED: gates {got}, wanted {want} every time")
            ok = False
        elif len(set(got)) != 1:
            print(f"   FAILED: nondeterministic, gates {got}")
            ok = False
        else:
            print(f"   ok: {args.k}/{args.k}")

    print("\n4. the shrink preserved the pathology")
    spec = json.loads((task / "verifier" / "cells.json").read_text())
    for fam, f in spec["families"].items():
        counts = [(n, c["cis_iterations"]) for n, c in f["cells"].items()]
        reached = [c for _, c in counts if c is not None]
        unreached = [n for n, c in counts if c is None]
        lo, hi = min(reached), max(reached)
        print(f"   {fam}: unaccelerated {lo} .. {hi} iterations where it converges"
              + (f", never within cap on {', '.join(unreached)}" if unreached else ""))
        if not unreached and hi / lo < 100:
            print(f"   FAILED: {fam} spread only {hi/lo:.0f}x and every cell converges -- "
                  "the shrink flattened it")
            ok = False
    if ok:
        print("   ok")

    print("\n5. wall-clock")
    worst = max(walls)
    print(f"   worst verification {worst:.0f}s of a {args.budget_seconds:.0f}s budget")
    if worst > args.budget_seconds:
        print("   FAILED: too slow for an agent to iterate against")
        ok = False

    print("\n" + ("TASK VALID" if ok else "TASK NOT VALID"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
