#!/usr/bin/env python3
"""Measure the experimental fixed-point repair ``scheme.defect_omega``.

The shipped GSIS blends the synthetic
solution ``M*`` with the kinetic moment ``M(f)``::

    M^{n+1} = beta M* + (1 - beta) M(f)

At a fixed point the bracket ``o = M* - M(f)`` is stationary and non-zero --
that offset *is* the CIS/GSIS discrepancy of docs/LIMITATIONS.md #1.

``defect_omega = w`` carries a relaxed running estimate of that offset and
subtracts it before the blend::

    M^{n+1} = M(f) + beta (o - d),    d <- (1-w) d + w o

At a fixed point ``d = o``, the blend vanishes and ``M = M(f)`` -- exactly
CIS's fixed-point condition, for any ``beta``.  ``w = 0`` recovers the
shipped scheme bit for bit.

The question this script answers is whether the acceleration survives the
change: the defect is lagged by one iteration, which can in principle cancel
the very speed-up it is bolted onto.

Usage::

    python tools/defect_study.py [--tau-r 1e-1 ...] [--out defect_study.json]
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


def _case(tau_r, deg, mesh, accflag, tol, tmax, omega=0.0, tau_thr=None, every=1, anderson=0):
    c = Case.from_yaml(CASE)
    c.scheme.accflag = accflag
    c.scheme.defect_omega = omega
    c.scheme.defect_every = every
    c.scheme.defect_anderson = anderson
    c.flow.tau_r = tau_r
    if tau_thr is not None:
        c.flow.tau_thr = tau_thr
    c.dg.deg = deg
    c.mesh.file = f"../meshes/{mesh}.msh"
    c.iteration.tol = tol
    c.iteration.tmax = tmax
    c.output.field = c.output.run_record = c.output.runtime_log = False
    return c


def run(tau_r, deg, mesh, accflag, tol, tmax, omega=0.0, tau_thr=None, every=1, anderson=0):
    t0 = time.perf_counter()
    s = Solver(_case(tau_r, deg, mesh, accflag, tol, tmax, omega, tau_thr, every, anderson))
    r = s.run()
    return {
        "iters": r.iterations,
        "converged": bool(r.converged),
        "true_res": float(s.ctx.true_residual(s.mom, s.vdf)),
        "mass": float(r.mass),
        "wall": time.perf_counter() - t0,
        "temp": r.temp,
        "defect_norm": r.diagnostics.get("defect_norm"),
        "defect_steps": (s.acc.defect_steps if s.acc is not None else []),
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tau-r", type=float, nargs="+", default=[1e-1])
    p.add_argument("--omega", type=float, nargs="+",
                   default=[0.0, 0.1, 0.3, 0.5, 1.0])
    p.add_argument("--deg", type=int, default=3)
    p.add_argument("--mesh", default="A1_Nx11_Ny11")
    p.add_argument("--tol", type=float, default=1e-13)
    p.add_argument("--tmax", type=int, default=200_000)
    p.add_argument("--every", type=int, default=1,
                   help="0 = freeze the defect until the inner iteration converges")
    p.add_argument("--anderson", type=int, default=0,
                   help="Anderson window on the outer defect sequence")
    p.add_argument("--out", default="defect_study.json")
    args = p.parse_args(argv)

    rows = []
    for tau_r in args.tau_r:
        cis = run(tau_r, args.deg, args.mesh, 0, args.tol, args.tmax)
        tc = cis["temp"]
        print(f"\ntau_R = {tau_r:.0e}   CIS: {cis['iters']} iters, "
              f"true_res {cis['true_res']:.2e}, mass {cis['mass']:.10f}")
        hdr = (f"{'omega':>6} {'GSIS it':>8} {'conv':>5} {'gap vs CIS':>11} "
               f"{'true_res':>10} {'mass':>14} {'|d|':>9} {'wall/s':>7}")
        print(hdr)
        print("-" * len(hdr))
        for w in args.omega:
            g = run(tau_r, args.deg, args.mesh, 1, args.tol, args.tmax,
                    omega=w, every=args.every, anderson=args.anderson)
            gap = float(np.abs(tc - g["temp"]).max() / np.abs(tc).max())
            print(f"{w:6.2f} {g['iters']:8d} {str(g['converged']):>5} "
                  f"{gap:11.3e} {g['true_res']:10.2e} {g['mass']:14.10f} "
                  f"{(g['defect_norm'] or 0.0):9.2e} {g['wall']:7.2f}"
                  f"   outer={len(g['defect_steps'])}", flush=True)
            rows.append({"tau_r": tau_r, "omega": w, "gap": gap,
                         "cis_iters": cis["iters"], "cis_mass": cis["mass"],
                         "cis_true_res": cis["true_res"],
                         **{k: v for k, v in g.items() if k != "temp"}})

    Path(args.out).write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
