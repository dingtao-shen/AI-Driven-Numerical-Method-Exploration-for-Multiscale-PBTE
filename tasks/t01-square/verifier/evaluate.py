"""Score one case in the *pristine* tree.

Two jobs, both done with code the submission never had access to:

1. **Correctness.**  The transport residual is computed from the submitted
   distribution here, not by the submitted solver -- a ``true_residual``
   that returns zero would otherwise pass everything.  A diffusely
   reflecting wall's emission is part of the operator and is rebuilt from
   the submitted distribution first.

2. **The unit of cost.**  The score is in *sweep-equivalents*: wall-clock
   divided by the time of one fine transport sweep, measured here, warm, in
   the same session as the submission's run.  The unaccelerated solver's
   cost is re-timed here too (its setup, and its seconds per iteration from
   a short forced run) so that both sides of the ratio are measured on the
   same machine under the same load; only its iteration *count* comes from
   calibration, since that is deterministic.
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np

pristine, case_path, dump, ref = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
sys.path.insert(0, pristine)

from pybte import Case, Solver                                    # noqa: E402
from pybte.moments import compute_moments                         # noqa: E402

d = np.load(dump)
r = np.load(ref)

case = Case.from_yaml(case_path)
case.output.field = case.output.run_record = case.output.runtime_log = False
case.iteration.tmax = 0                       # build the operators, iterate zero times
solver = Solver(case)

vdf = d["vdf"]
out = {}
if vdf.shape != solver.vdf.shape:
    out["error"] = (f"vdf shape {vdf.shape} does not match the discretisation "
                    f"{solver.vdf.shape}")
    print(json.dumps(out))
    raise SystemExit(0)

# -- correctness ---------------------------------------------------------------
solver.vdf[:] = vdf
compute_moments(solver.vdf, solver.ctx.cxv, solver.ctx.cyv, solver.ctx.domega,
                solver.integrals.int_tri, case.flow.cv, solver.mom)
if solver.bcdata.has_nonthermalising:
    solver.bcdata.update_wall_flux(solver.vdf, solver.ctx.cxv, solver.ctx.cyv,
                                   solver.ctx.domega, solver.integrals.int_tri_fc,
                                   solver.mesh)
out["transport_residual"] = float(solver.ctx.true_residual(solver.mom, solver.vdf))

ref_t = r["temp"]
got_t = d["temp"]
if got_t.shape != ref_t.shape:
    out["error"] = f"temp shape {got_t.shape} != reference {ref_t.shape}"
    print(json.dumps(out))
    raise SystemExit(0)
scale = float(np.abs(ref_t).max())
out["field_gap"] = float(np.abs(got_t - ref_t).max() / scale)
out["field_gap_from_vdf"] = float(np.abs(solver.mom.temp - ref_t).max() / scale)
out["reference_transport_residual"] = float(r["transport_residual"])

# -- the unit: one warm fine sweep, median of five ---------------------------------
solver.ctx.sweep(solver.mom, solver.vdf)                      # warm
ts = []
for _ in range(5):
    t0 = time.perf_counter()
    solver.ctx.sweep(solver.mom, solver.vdf)
    ts.append(time.perf_counter() - t0)
out["t_sweep_s"] = float(np.median(ts))

# -- the unaccelerated solver's cost, re-timed now ---------------------------------
c2 = Case.from_yaml(case_path)
c2.output.field = c2.output.run_record = c2.output.runtime_log = False
c2.iteration.tol = 0.0                         # never satisfied: run exactly tmax
c2.iteration.tmax = 200
t0 = time.perf_counter()
s2 = Solver(c2)
out["cis_setup_s"] = time.perf_counter() - t0
t0 = time.perf_counter()
rec2 = s2.run()
out["cis_s_per_iter"] = (time.perf_counter() - t0) / max(1, rec2.iterations)
print(json.dumps(out))
