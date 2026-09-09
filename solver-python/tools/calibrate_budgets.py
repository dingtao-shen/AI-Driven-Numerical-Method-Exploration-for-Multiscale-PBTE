#!/usr/bin/env python3
"""Calibrate the per-case iteration budgets for the benchmark ladder.

A task budget has to satisfy two things at once:

* the fixed-point-preserving oracle must clear it with headroom;
* the unaccelerated solver must not come close.

Both numbers have to be measured on the *shrunk* case the task actually
ships, not on the research discretisation -- a shrink that accidentally lets
CIS converge quickly destroys the task, so this script also checks that the
pathology survives shrinking.

Three schemes per Knudsen number:

``cis``       ``accflag=0``.  Correct but slow; the thing the task removes.
``gsis``      ``accflag=1`` as shipped.  Fast, but converges to a different
              discrete fixed point, so it cannot be the oracle.
``oracle``    ``accflag=1`` plus the defect correction with Anderson on the
              outer loop.  Preserves CIS's fixed point.

Usage::

    python tools/calibrate_budgets.py [--deg 2] [--npole 10] [--nazim 20]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from pybte import Case, Solver                                    # noqa: E402

CASE = HERE.parent / "cases" / "cavity_tauR1e-3_cis.yaml"

SCHEMES = {
    "cis":    dict(accflag=0, omega=0.0, anderson=0),
    "gsis":   dict(accflag=1, omega=0.0, anderson=0),
    "oracle": dict(accflag=1, omega=1.0, anderson=8),
}


def build(tau_r, scheme, args, tmax):
    c = Case.from_yaml(CASE)
    s = SCHEMES[scheme]
    c.scheme.accflag = s["accflag"]
    c.scheme.defect_omega = s["omega"]
    c.scheme.defect_anderson = s["anderson"]
    c.scheme.defect_every = 0
    c.flow.tau_r = tau_r
    c.flow.tau_n = args.tau_n
    c.dg.deg = args.deg
    c.velmesh.npole = args.npole
    c.velmesh.nazim = args.nazim
    c.mesh.file = f"../meshes/{args.mesh}.msh"
    c.iteration.tol = args.tol
    c.iteration.tmax = tmax
    c.output.field = c.output.run_record = c.output.runtime_log = False
    return c


def run(tau_r, scheme, args, tmax):
    t0 = time.perf_counter()
    s = Solver(build(tau_r, scheme, args, tmax))
    r = s.run()
    return {
        "scheme": scheme, "tau_r": tau_r,
        "iterations": r.iterations, "converged": bool(r.converged),
        "mass": float(r.mass),
        "true_residual": float(s.ctx.true_residual(s.mom, s.vdf)),
        "outer": len(s.acc.defect_steps) if s.acc is not None else 0,
        "wall": time.perf_counter() - t0,
        "temp": r.temp,
        "n_tris": s.n_tris, "ndir": int(s.ctx.cxv.size),
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mesh", default="A1_Nx11_Ny11")
    p.add_argument("--deg", type=int, default=2)
    p.add_argument("--npole", type=int, default=10)
    p.add_argument("--nazim", type=int, default=20)
    p.add_argument("--tau-n", type=float, default=1e5)
    p.add_argument("--tol", type=float, default=1e-8)
    p.add_argument("--tau-r", type=float, nargs="+",
                   default=[1e0, 1e-1, 1e-2, 1e-3, 1e-4])
    p.add_argument("--cis-tmax", type=int, default=200_000)
    p.add_argument("--acc-tmax", type=int, default=100_000)
    p.add_argument("--out", default="docs/budgets.json")
    args = p.parse_args(argv)

    rows = []
    hdr = (f"{'tau_R':>7} {'scheme':>7} {'iters':>8} {'conv':>5} {'outer':>6} "
           f"{'true_res':>10} {'mass':>14} {'gap vs CIS':>11} {'wall/s':>7}")
    print(f"# shrunk case: {args.mesh}, DEG={args.deg}, "
          f"{args.npole}x{args.nazim} angles, tol={args.tol:.0e}")
    print(hdr); print("-" * len(hdr), flush=True)

    for tau_r in args.tau_r:
        ref = None
        for scheme in ("cis", "gsis", "oracle"):
            tmax = args.cis_tmax if scheme == "cis" else args.acc_tmax
            r = run(tau_r, scheme, args, tmax)
            if scheme == "cis":
                ref = r["temp"] if r["converged"] else None
            gap = (float(np.abs(ref - r["temp"]).max() / np.abs(ref).max())
                   if ref is not None else float("nan"))
            r["gap_vs_cis"] = gap
            print(f"{tau_r:7.0e} {scheme:>7} {r['iterations']:8d} "
                  f"{str(r['converged']):>5} {r['outer']:6d} "
                  f"{r['true_residual']:10.2e} {r['mass']:14.10f} "
                  f"{gap:11.3e} {r['wall']:7.1f}", flush=True)
            rows.append({k: v for k, v in r.items() if k != "temp"})

    Path(args.out).write_text(json.dumps(
        {"config": vars(args), "rows": rows}, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
