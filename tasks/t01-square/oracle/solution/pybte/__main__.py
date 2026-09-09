"""Command-line interface.

    python -m pybte run     case.yaml [--tmax N] [--tol T] [--quiet]
    python -m pybte convert control.in [-o case.yaml] [--mesh-dir DIR]
    python -m pybte info    case.yaml
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def _cmd_run(args) -> int:
    from .config import Case
    from .driver import Solver

    case = Case.from_yaml(args.case)
    if args.tmax is not None:
        case.iteration.tmax = args.tmax
    if args.tol is not None:
        case.iteration.tol = args.tol
    if args.true_residual:
        case.iteration.true_residual = True
    if args.out is not None:
        case.output.dir = args.out
    if args.variant is not None:
        case.scheme.acc_variant = args.variant
    if args.no_output:
        case.output.field = case.output.run_record = case.output.runtime_log = False

    solver = Solver(case)
    if not args.quiet:
        print(f"pybte: {case.mesh_path().name}  N_TRIS={solver.n_tris} "
              f"N_FCS={solver.mesh.n_faces}  DEG={case.dg.deg}  "
              f"NPOLE={case.velmesh.npole} NAZIM={case.velmesh.nazim}")
        print(f"       scheme={'GSIS ' + case.scheme.acc_variant if case.scheme.accflag else 'CIS'}"
              f"  tau_R={case.flow.tau_r:g} tau_N={case.flow.tau_n:g} "
              f"tau_C={case.tau_c:g}  tol={case.iteration.tol:g}")
        print(f"       setup {solver.setup_seconds:.2f}s, sweep mode {solver.ctx.mode}")

    rec = solver.run(progress_every=args.progress)
    written = solver.write_outputs(rec)

    if not args.quiet:
        print(f"  iterations      {rec.iterations}"
              f"{'  (converged)' if rec.converged else '  (TMAX reached)'}")
        print(f"  residual        {rec.final_residual:.6e}")
        print(f"  mass (int T dA) {rec.mass:.12f}")
        print(f"  sweeps          {rec.sweep_count}")
        print(f"  factorisations  {rec.factorisation_count}")
        print(f"  wall clock      {rec.wall_clock:.3f} s"
              f"  ({rec.wall_clock / max(rec.iterations, 1) * 1e3:.2f} ms/iter)")
        for k, v in written.items():
            print(f"  wrote {k:12s} {v}")
    return 0 if rec.converged else 1


def _cmd_convert(args) -> int:
    from .config import case_from_control_in

    case = case_from_control_in(args.control, mesh_dir=args.mesh_dir)
    out = Path(args.out) if args.out else Path(args.control).with_suffix(".yaml")
    case.to_yaml(out)
    print(f"wrote {out}")
    return 0


def _cmd_info(args) -> int:
    from .config import Case
    from .driver import Solver

    case = Case.from_yaml(args.case)
    s = Solver(case)
    print(f"mesh            {case.mesh_path()}")
    print(f"  nodes         {s.mesh.n_nodes}")
    print(f"  triangles     {s.mesh.n_tris}")
    print(f"  faces         {s.mesh.n_faces} ({s.mesh.n_faces_b} boundary, "
          f"{s.mesh.n_faces_i} interior)")
    print(f"  hmin          {s.mesh.hmin:.6e}")
    print(f"  area          {s.mesh.tri_area.sum():.12f}")
    print(f"discretisation  DEG={case.dg.deg} NDOF_TRI={case.ndof_tri} "
          f"NDOF_FC={case.ndof_fc}")
    print(f"  directions    {s.ndir} (sum dOmega = {s.vel.domega.sum():.12f}, "
          f"4pi = {4 * np.pi:.12f})")
    print(f"  operators     {s.ctx.mode}, "
          f"{s.ctx.operator_bytes(s.n_tris, case.ndof_tri) / 1e6:.0f} MB")
    print("boundaries")
    for b in case.boundaries:
        n = int(np.count_nonzero(s.mesh.face_bc == case.boundaries.index(b)))
        print(f"  {b.name:8s} phyid={b.phyid:3d} type={b.type:16s} "
              f"T={b.temp:g}  faces={n}")
    if s.acc is not None:
        print(f"GSIS variant {case.scheme.acc_variant}: "
              f"global matrix {s.acc.K.shape[0]} x {s.acc.K.shape[0]}, "
              f"nnz {s.acc.K.nnz}")
    print(f"config hash     {case.config_hash}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="pybte", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="solve a case")
    r.add_argument("case")
    r.add_argument("--tmax", type=int, default=None)
    r.add_argument("--tol", type=float, default=None)
    r.add_argument("--out", default=None, help="override output.dir")
    r.add_argument("--variant", default=None, choices=["A", "B"])
    r.add_argument("--true-residual", action="store_true")
    r.add_argument("--progress", type=int, default=0,
                   help="print every N iterations (0 = silent)")
    r.add_argument("--no-output", action="store_true")
    r.add_argument("--quiet", action="store_true")
    r.set_defaults(fn=_cmd_run)

    c = sub.add_parser("convert", help="control.in -> case.yaml")
    c.add_argument("control")
    c.add_argument("-o", "--out", default=None)
    c.add_argument("--mesh-dir", default=None)
    c.set_defaults(fn=_cmd_convert)

    i = sub.add_parser("info", help="describe a case without solving it")
    i.add_argument("case")
    i.set_defaults(fn=_cmd_info)

    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
