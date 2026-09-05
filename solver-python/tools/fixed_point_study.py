#!/usr/bin/env python3
"""Quantify the CIS/GSIS fixed-point discrepancy.

Proposal 1 §10 asserts, as "the load-bearing property for the whole
benchmark", that CIS and GSIS converge to the same *discrete* fixed point to
``rtol=1e-8``.  They do not, in the reference or in this port.  This script
measures by how much, and how that gap scales, so that Proposal 2 can be
written against what the solver actually does.

Three numbers per configuration:

``true_res``    ``||A f - b|| / ||b||`` of the discrete transport system,
                evaluated with no sweep.  Zero exactly at the kinetic fixed
                point.  CIS drives it to round-off; GSIS does not.
``gap``         ``||T_CIS - T_GSIS||_inf / ||T_CIS||_inf`` after both have
                converged on the iterate residual.
``consistency`` ``||UQ_T - T_VDF||_inf / ||T_VDF||_inf`` for one macroscopic
                solve applied to the *converged CIS* solution.  This is the
                defect at its source: if it were zero, the CIS fixed point
                would be a GSIS fixed point.

Usage::

    python tools/fixed_point_study.py [--out fixed_point_study.json]
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
from pybte.acceleration.hot_source import hot_source              # noqa: E402
from pybte.acceleration.local_solve import local_solve            # noqa: E402

CASE = HERE.parent / "cases" / "cavity_tauR1e-3_cis.yaml"


def _case(tau_r, deg, mesh, accflag, tol, tmax):
    c = Case.from_yaml(CASE)
    c.scheme.accflag = accflag
    c.flow.tau_r = tau_r
    c.dg.deg = deg
    c.mesh.file = f"../meshes/{mesh}.msh"
    c.iteration.tol = tol
    c.iteration.tmax = tmax
    c.output.field = c.output.run_record = c.output.runtime_log = False
    return c


def consistency_defect(cis_solver, tau_r, deg, mesh):
    """One macroscopic solve applied to the converged kinetic solution."""
    cg = _case(tau_r, deg, mesh, 1, 1e-13, 1)
    sg = Solver(cg)
    sg.vdf[:] = cis_solver.vdf
    acc = sg.acc
    acc.aa_src, acc._hot_work = hot_source(
        sg.vdf, acc.gx_t, acc.gy_t, sg.ctx.cxv, sg.ctx.cyv, sg.ctx.domega,
        cg.tau_c, acc.nd, acc.inv_aa_sol, acc.scale, acc._hot_work)
    ffa = acc.gsolver.rhs(acc.aa_src, sg.vdf, sg.ctx.cxv, sg.ctx.cyv, sg.ctx.domega)
    u_trace = acc.gsolver.solve(ffa)
    uq = local_solve(acc.aa_src, acc.aa_trace, u_trace, sg.mesh.tri_faces)

    nd = acc.nd
    t_acc = uq[0:nd, :]
    t_vdf = cis_solver.mom.ts
    d = np.abs(t_acc - t_vdf) / np.abs(t_vdf).max()

    cen = sg.mesh.nodes[sg.mesh.tri_nodes].mean(axis=1)
    dist = np.minimum.reduce([cen[:, 0], 1 - cen[:, 0], cen[:, 1], 1 - cen[:, 1]])
    per_el = d.max(axis=0)
    return {
        "consistency_inf": float(per_el.max()),
        "consistency_bulk": float(np.median(per_el[dist > 0.2])) if (dist > 0.2).any() else None,
        "consistency_wall": float(np.median(per_el[dist < 0.02])) if (dist < 0.02).any() else None,
    }


def one(tau_r, deg, mesh, tol=1e-13, tmax=200_000):
    row = {"tau_r": tau_r, "deg": deg, "mesh": mesh}
    solvers = {}
    for accflag, key in ((0, "cis"), (1, "gsis")):
        t0 = time.perf_counter()
        s = Solver(_case(tau_r, deg, mesh, accflag, tol, tmax))
        r = s.run()
        row[f"{key}_iters"] = r.iterations
        row[f"{key}_converged"] = bool(r.converged)
        row[f"{key}_true_res"] = float(s.ctx.true_residual(s.mom, s.vdf))
        row[f"{key}_mass"] = float(r.mass)
        row[f"{key}_wall"] = time.perf_counter() - t0
        solvers[key] = (s, r)
    tc = solvers["cis"][1].temp
    tg = solvers["gsis"][1].temp
    row["gap"] = float(np.abs(tc - tg).max() / np.abs(tc).max())
    row["n_tris"] = solvers["cis"][0].n_tris
    row.update(consistency_defect(solvers["cis"][0], tau_r, deg, mesh))
    return row


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="fixed_point_study.json")
    args = p.parse_args(argv)

    configs = [
        # Knudsen sweep, fixed discretisation
        (1.0, 3, "A1_Nx11_Ny11"),
        (1e-1, 3, "A1_Nx11_Ny11"),
        (1e-2, 3, "A1_Nx11_Ny11"),
        # polynomial order, fixed Kn and mesh
        (1e-1, 1, "A1_Nx11_Ny11"),
        (1e-1, 2, "A1_Nx11_Ny11"),
        # mesh refinement, fixed Kn and order
        (1e-1, 3, "A1_Nx6_Ny6"),
        (1e-1, 3, "A1_Nx21_Ny21"),
    ]
    rows = []
    hdr = (f"{'tau_R':>7} {'DEG':>3} {'N_TRIS':>7} | {'CIS it':>7} {'CIS true':>10} | "
           f"{'GSIS it':>7} {'GSIS true':>10} | {'gap':>9} {'consist':>9} {'bulk':>9}")
    print(hdr)
    print("-" * len(hdr))
    for tau_r, deg, mesh in configs:
        row = one(tau_r, deg, mesh)
        rows.append(row)
        print(f"{row['tau_r']:7.0e} {row['deg']:3d} {row['n_tris']:7d} | "
              f"{row['cis_iters']:7d} {row['cis_true_res']:10.2e} | "
              f"{row['gsis_iters']:7d} {row['gsis_true_res']:10.2e} | "
              f"{row['gap']:9.2e} {row['consistency_inf']:9.2e} "
              f"{row['consistency_bulk']:9.2e}", flush=True)

    Path(args.out).write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
