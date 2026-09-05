#!/usr/bin/env python3
"""Compare a pybte run against a Fortran stage dump.

This is the primary debugging instrument of the port: it localises any
discrepancy to a single subroutine instead of leaving you to bisect an
iteration count.

    python tools/compare_stages.py DUMPDIR [--case CASE.yaml] [--iters 1,2,3,10]

Without ``--case`` the case is reconstructed from ``DUMPDIR/../control.in``,
so a dump directory is self-describing.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

from dump_fortran_stages import FortranDump          # noqa: E402
from pybte import Case, Solver                        # noqa: E402
from pybte.config import case_from_control_in         # noqa: E402


class Report:
    def __init__(self, rtol_default=1e-12):
        self.rows = []
        self.rtol_default = rtol_default
        self.failures = 0

    def check(self, name, got, want, rtol=None, exact=False):
        rtol = self.rtol_default if rtol is None else rtol
        got = np.asarray(got)
        want = np.asarray(want)
        if got.shape != want.shape:
            self.rows.append((name, "SHAPE", f"{got.shape} vs {want.shape}", ""))
            self.failures += 1
            return False
        if exact:
            ok = np.array_equal(got, want)
            self.rows.append((name, "OK" if ok else "FAIL", "exact",
                              "bit-exact" if ok else ""))
            self.failures += 0 if ok else 1
            return ok
        err = np.abs(got - want)
        scale = np.abs(want).max()
        maxabs = float(err.max())
        rel = maxabs / scale if scale else maxabs
        ok = rel <= rtol
        bit = "bit-exact" if np.array_equal(got, want) else ""
        self.rows.append((name, "OK" if ok else "FAIL", f"rel={rel:.3e}", bit))
        self.failures += 0 if ok else 1
        return ok

    def print(self):
        w = max((len(r[0]) for r in self.rows), default=10)
        for name, status, detail, extra in self.rows:
            print(f"  {name:<{w}}  {status:4s}  {detail:20s} {extra}")
        print(f"  --> {len(self.rows) - self.failures}/{len(self.rows)} checks passed")


def compare(dumpdir, case=None, iters=(1, 2, 3, 10), rtol=1e-12, quiet=False,
            variant=None):
    d = FortranDump(dumpdir)
    if case is None:
        control = Path(dumpdir).parent / "control.in"
        case = case_from_control_in(control, mesh_dir=str(Path(dumpdir).parent))
    if variant:
        case.scheme.acc_variant = variant
    case.iteration.tmax = max(iters)
    case.restart.error_on_stale = False

    rep = Report(rtol)
    solver = Solver(case)

    # -- setup ---------------------------------------------------------
    rep.check("NODES", solver.mesh.nodes, d["NODES"])
    rep.check("TRIANGLES_TAG", solver.mesh.fortran_triangles_tag(),
              d["TRIANGLES_TAG"], exact=True)
    rep.check("TRIANGLES_INF", solver.mesh.fortran_triangles_inf(), d["TRIANGLES_INF"])
    rep.check("TRIANGLES_Hmin", solver.mesh.tri_hmin, d["TRIANGLES_Hmin"])
    rep.check("FACES_TAG", solver.mesh.fortran_faces_tag(), d["FACES_TAG"], exact=True)
    rep.check("FACES_INF", solver.mesh.face_len, d["FACES_INF"])
    rep.check("THE", solver.vel.the, d["THE"])
    rep.check("PHI", solver.vel.phi, d["PHI"])
    rep.check("DOMEGA", solver.vel.domega, d["DOMEGA"])
    rep.check("CX", solver.vel.cx, d["CX"])
    rep.check("CY", solver.vel.cy, d["CY"])
    rep.check("NODFUN_TRI_REF", solver.basis.nodfun_tri, d["NODFUN_TRI_REF"])
    rep.check("NODFUN_FC_REF", solver.basis.nodfun_fc, d["NODFUN_FC_REF"])
    rep.check("TRI_ORDER", solver.order.fortran_tri_order(), d["TRI_ORDER"], exact=True)

    I = solver.integrals
    rep.check("INT_NODFUNC_TRI", I.int_tri, d["INT_NODFUNC_TRI"])
    rep.check("INT_NODFUNC_TRI_TRI", I.int_tri_tri, d["INT_NODFUNC_TRI_TRI"])
    rep.check("INT_NODFUNC_TRI_TRI_X", I.int_tri_tri_x, d["INT_NODFUNC_TRI_TRI_X"])
    rep.check("INT_NODFUNC_TRI_TRI_Y", I.int_tri_tri_y, d["INT_NODFUNC_TRI_TRI_Y"])
    rep.check("INT_NODFUNC_TRI_FC", I.int_tri_fc, d["INT_NODFUNC_TRI_FC"])
    rep.check("INT_NODFUNC_TRI_FC_FC", I.int_tri_fc_fc, d["INT_NODFUNC_TRI_FC_FC"])
    rep.check("INT_NODFUNC_FC_FC", I.int_fc_fc, d["INT_NODFUNC_FC_FC"])
    rep.check("INT_NODFUNC_FC", I.int_fc, d["INT_NODFUNC_FC"])
    rep.check("INT_NODFUNC_TRI_TRI_FC", I.int_tri_tri_fc, d["INT_NODFUNC_TRI_TRI_FC"])
    rep.check("INT_NODFUNC_TRI_TRI_TRI", I.int_tri_tri_tri(solver.mesh.tri_area),
              d["INT_NODFUNC_TRI_TRI_TRI"])
    rep.check("NODFUN_QUA_P", I.nodfun_qua_p, d["NODFUN_QUA_P"])

    gsis = solver.acc is not None
    if gsis:
        acc = solver.acc
        rep.check("ST", acc.st, d["ST"])
        rep.check("inv_AA_SOL", acc.inv_aa_sol_fortran(), d["inv_AA_SOL"])
        rep.check("AA_TRACE", acc.aa_trace_fortran(), d["AA_TRACE"])
        rep.check("BA_SOL", acc.ba_sol_fortran(), d["BA_SOL"])
        K = d.csr()
        rep.check("K (global matrix)", acc.global_matrix().toarray(), K.toarray())

    # -- per-iteration stages -------------------------------------------
    for step in range(1, max(iters) + 1):
        solver.step()
        if step not in iters:
            continue
        if d.has_iter(step, "vdf_after_sweep") and not gsis:
            rep.check(f"it{step}:VDF", solver.vdf_fortran(),
                      d.iter_array(step, "vdf_after_sweep"))
        if gsis:
            rep.check(f"it{step}:AA_SRC", solver.acc.aa_src, d.iter_array(step, "AA_SRC"))
            rep.check(f"it{step}:FFA", solver.acc.ffa.ravel(order="F"), d.iter_array(step, "FFA"))
            rep.check(f"it{step}:U_TRACE", solver.acc.u_trace, d.iter_array(step, "U_TRACE"))
            rep.check(f"it{step}:UQ", solver.acc.uq, d.iter_array(step, "UQ"))
            rep.check(f"it{step}:VDF_corr", solver.vdf_fortran(),
                      d.iter_array(step, "vdf_after_correction"))
            tag = "after_correction"
        else:
            tag = "after_moments"
        rep.check(f"it{step}:T_s", solver.mom.ts, d.iter_array(step, f"Ts_{tag}"))
        rep.check(f"it{step}:Qx_s", solver.mom.qxs, d.iter_array(step, f"Qxs_{tag}"))
        rep.check(f"it{step}:Qy_s", solver.mom.qys, d.iter_array(step, f"Qys_{tag}"))
        rep.check(f"it{step}:Temp", solver.mom.temp, d.iter_array(step, f"Temp_{tag}"))

    if not quiet:
        rep.print()
    return rep


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("dumpdir")
    p.add_argument("--case", default=None)
    p.add_argument("--iters", default="1,2,3,10")
    p.add_argument("--rtol", type=float, default=1e-12)
    p.add_argument("--variant", default=None, choices=["A", "B"])
    args = p.parse_args(argv)

    case = Case.from_yaml(args.case) if args.case else None
    iters = tuple(int(v) for v in args.iters.split(","))
    rep = compare(args.dumpdir, case, iters, args.rtol, variant=args.variant)
    return 1 if rep.failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
