"""Reproduce the reference publication's 2-D square table at its own discretisation.

Five (Kn_R, Kn_N) points, 200 triangles, P3, the published angular
resolutions, the published stopping rule (L2-relative iterate residual
< 1e-7).  Both the unaccelerated and the accelerated iteration counts are
compared with the published ones, and the field gap between the two
schemes' converged solutions is reported: it is the fixed-point
displacement of docs/LIMITATIONS.md #1, measured on the published test
points themselves.

Usage::

    python tools/reproduce_published_table.py
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pybte import Case, Solver
ROWS = [  # (Kn_R, Kn_N, npole, nazim, paper_GSIS, paper_CIS)
    (1e-3, 1e5, 20, 40, 43, None),
    (1e-2, 1e5, 40, 80, 24, 13234),
    (1e-1, 1e5, 40, 80, 26, 269),
    (1e0,  1e0, 80, 160, 28, 29),
    (1e1,  1e-2, 40, 80, 47, 1883),
]
def mk(kr, kn, npole, nazim, accflag, tmax):
    c = Case.from_yaml(Path(__file__).resolve().parent.parent / "cases" / "cavity_tauR1e-3_cis.yaml")
    c.scheme.accflag = accflag; c.flow.tau_r = kr; c.flow.tau_n = kn
    c.flow.tau_thr = 100.0 if kr == 10 else 1.0
    c.dg.deg = 3; c.velmesh.npole = npole; c.velmesh.nazim = nazim
    c.mesh.file = "../meshes/A1_Nx11_Ny11.msh"
    c.iteration.tol = 1e-7; c.iteration.tmax = tmax
    c.output.field = c.output.run_record = c.output.runtime_log = False
    return c
print(f"{'Kn_R':>6} {'Kn_N':>6} {'angles':>7} | {'GSIS n':>7} {'paper':>6} {'wall/s':>7} {'tr':>9} | "
      f"{'CIS n':>7} {'paper':>6} {'wall/s':>7} {'tr':>9} | {'gap':>9}")
print("-"*110)
for kr, kn, npole, nazim, pg, pc in ROWS:
    t0 = time.perf_counter(); sg = Solver(mk(kr, kn, npole, nazim, 1, 2000)); rg = sg.run(); wg = time.perf_counter()-t0
    trg = sg.ctx.true_residual(sg.mom, sg.vdf)
    line = (f"{kr:6.0e} {kn:6.0e} {npole:3d}x{nazim:<3d} | {rg.iterations:6d}{'' if rg.converged else '*'} "
            f"{pg:6d} {wg:7.0f} {trg:9.2e} | ")
    if pc is None:
        line += f"{'skip':>7} {'>1M':>6} {'':>7} {'':>9} | {'':>9}"
    else:
        t0 = time.perf_counter(); sc = Solver(mk(kr, kn, npole, nazim, 0, 40000)); rc = sc.run(); wc = time.perf_counter()-t0
        trc = sc.ctx.true_residual(sc.mom, sc.vdf)
        import numpy as np
        gap = float(np.abs(rg.temp-rc.temp).max()/np.abs(rc.temp).max())
        line += (f"{rc.iterations:6d}{'' if rc.converged else '*'} {pc:6d} {wc:7.0f} {trc:9.2e} | {gap:9.2e}")
    print(line, flush=True)
print("\n* = not converged within tmax.  tol = 1e-7 L2-relative on the iterate, as in the paper's Eq (52).")
