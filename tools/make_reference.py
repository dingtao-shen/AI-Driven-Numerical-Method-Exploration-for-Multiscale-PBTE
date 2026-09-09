#!/usr/bin/env python3
"""Generate a task's reference solutions, certified by transport residual.

**A reference is certified by ``||A f - b|| / ||b||``, never by an iterate
tolerance.**  The iterate residual measures the step between successive
iterates, not the distance to the answer, and in the diffusive regime the two
differ by orders of magnitude.  Worse, the fixed-point-preserving scheme's
outer loop has a round-off floor of its own, so at some Knudsen numbers it
cannot meet a tight *iterate* tolerance while sitting comfortably at a tight
*transport* residual.  Certifying on the iterate residual would silently ship
references worse than the runs they grade.

Produced by the full solver at the task's exact discretisation, so the
reference and the environment differ only in what the ablation removed.

Usage::

    python tools/make_reference.py tasks/t01-diffusive-acceleration
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", type=Path)
    ap.add_argument("--tol", type=float, default=1e-10,
                    help="iterate tolerance for the reference run -- two decades "
                         "tighter than the runs it grades")
    ap.add_argument("--certify", type=float, default=1e-8,
                    help="required transport residual; an order of magnitude "
                         "below the verifier's own threshold")
    ap.add_argument("--tmax", type=int, default=60_000)
    ap.add_argument("--cis-tmax", type=int, default=40_000,
                    help="budget for the independent cross-check path")
    ap.add_argument("--cross-check", type=float, default=1e-5,
                    help="required agreement between the two code paths")
    args = ap.parse_args(argv)

    task = args.task.resolve()
    spec = yaml.safe_load((task / "ablation.yaml").read_text())
    src = (task / spec["source"]).resolve()
    sys.path.insert(0, str(src))
    from pybte import Case, Solver                                # noqa: E402

    env_cases = sorted((task / "environment" / "cases").glob("*.yaml"))
    out = task / "reference"
    out.mkdir(exist_ok=True)
    index = {}

    print(f"{'case':>10} {'iters':>7} {'outer':>6} {'transport res':>14} "
          f"{'int T dA':>14} {'wall/s':>7}  certified")
    for cf in env_cases:
        c = Case.from_yaml(cf)
        # the reference is produced by the Krylov solve of (I - T) u = g:
        # the same discrete fixed point as source iteration, reached to
        # round-off rather than to a tolerance
        c.scheme.accflag = 0
        c.scheme.method = "krylov"
        c.iteration.tol = args.tol
        c.iteration.tmax = args.tmax
        c.output.field = c.output.run_record = c.output.runtime_log = False
        c.mesh.file = str((cf.parent / c.mesh.file).resolve())

        t0 = time.perf_counter()
        s = Solver(c)
        r = s.run()
        tr = float(s.ctx.true_residual(s.mom, s.vdf))
        wall = time.perf_counter() - t0
        ok = tr < args.certify
        print(f"{cf.stem:>10} {r.iterations:7d} {'-':>6} "
              f"{tr:14.3e} {r.mass:14.10f} {wall:7.1f}  {'yes' if ok else 'NO'}",
              flush=True)
        if not ok:
            print(f"\nrefusing to ship {cf.stem}: transport residual {tr:.3e} "
                  f"is above the {args.certify:.0e} certificate", file=sys.stderr)
            return 1

        np.savez_compressed(
            out / f"{cf.stem}.npz",
            temp=r.temp, qx=r.qx, qy=r.qy,
            temp_dofs=s.mom.ts, iterations=r.iterations,
            transport_residual=tr, mass=r.mass,
            config_hash=c.config_hash)
        entry = {"iterations": int(r.iterations), "transport_residual": tr,
                 "mass": float(r.mass), "config_hash": c.config_hash}

        # Cross-check against the unaccelerated solver wherever it can reach
        # the answer.  This is what licenses using the reference as the
        # fixed-point standard: two independent code paths, one of them the
        # baseline the task ships, agreeing to well inside the check
        # threshold.  Never ship a reference produced by one path alone.
        cc = Case.from_yaml(cf)
        cc.scheme.accflag = 0
        cc.iteration.tol = args.tol
        cc.iteration.tmax = args.cis_tmax
        cc.output.field = cc.output.run_record = cc.output.runtime_log = False
        cc.mesh.file = str((cf.parent / cc.mesh.file).resolve())
        sc = Solver(cc)
        rc = sc.run()
        if rc.converged:
            gap = float(np.abs(rc.temp - r.temp).max() / np.abs(rc.temp).max())
            entry["cis_iterations"] = int(rc.iterations)
            entry["cis_gap"] = gap
            entry["cis_transport_residual"] = float(
                sc.ctx.true_residual(sc.mom, sc.vdf))
            mark = "ok" if gap < args.cross_check else "FAILED"
            print(f"{'':>10} cross-check vs unaccelerated: {rc.iterations} iters, "
                  f"gap {gap:.3e}  {mark}", flush=True)
            if gap >= args.cross_check:
                print(f"\nrefusing to ship {cf.stem}: the two code paths differ "
                      f"by {gap:.3e}", file=sys.stderr)
                return 1
        else:
            entry["cis_iterations"] = None
            print(f"{'':>10} cross-check vs unaccelerated: not reachable in "
                  f"{args.cis_tmax} iterations", flush=True)
        index[cf.stem] = entry

    (out / "index.json").write_text(json.dumps(index, indent=2))
    print(f"\nwrote {len(index)} references to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
