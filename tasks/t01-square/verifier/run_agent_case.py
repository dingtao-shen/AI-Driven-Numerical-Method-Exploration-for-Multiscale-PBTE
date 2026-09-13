"""Run one case in the *agent's* tree and dump what the verifier needs.

Touches only the public interface the task promises -- ``Case.from_yaml``,
``Solver(case).run()``, ``solver.vdf`` -- so a legitimate rewrite of the
internals cannot break it.

Timing is the whole of ``Case.from_yaml + Solver(case) + run()`` for the
graded case, after a warm-up solve of the *same kind of case on a coarser
mesh*: the warm-up gets every kernel the graded run uses compiled and its
cache loaded, so the measurement is of the method and not of the JIT, while
a different mesh means nothing about the graded case can have been
precomputed or memoised.
Setup is inside the timing on purpose -- a method that moves all its work
into "setup" still pays for it.

The warm-up must succeed on every family.  It is timed out of the score, so
a family whose warm-up throws is graded cold and charged for the JIT; the
result is reported as ``warm_up`` and ``checks.py`` refuses to score a cell
whose warm-up failed.
"""
from __future__ import annotations

import copy
import json
import resource
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import yaml

env, case_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, env)

from pybte import Case, Solver                                    # noqa: E402

spec = yaml.safe_load(Path(case_path).read_text())

# -- warm-up: the graded case's own wall types, on the coarse mesh of the ----
# -- same family, three iterations ------------------------------------------
# Every kernel the graded run will touch -- thermalising, reflecting or
# periodic walls -- is compiled and its disk cache loaded here, where it is
# not timed.  A warm-up that only exercised thermalising walls left the first
# reflecting-wall cell of a run to compile inside its own timing (~1 s, a
# fivefold penalty on a 0.2 s run), and one that dropped the periodic offsets
# while keeping the Master/Slave names failed the mesh reader's pairing and
# graded the whole periodic family cold.
warm_note = ""
try:
    w = copy.deepcopy(spec)
    mesh = Path(spec["mesh"]["file"])
    w["mesh"]["file"] = str(mesh.parent / f"{mesh.name.split('_')[0]}_Nx6_Ny6.msh")
    w["flow"]["tau_r"] = 1.0
    w["iteration"]["tmax"] = 3
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as fh:
        yaml.safe_dump(w, fh, sort_keys=False)
        wpath = fh.name
    wc = Case.from_yaml(wpath)
    wc.output.field = wc.output.run_record = wc.output.runtime_log = False
    Solver(wc).run()
except Exception as e:                                            # noqa: BLE001
    warm_note = f"warm-up failed: {type(e).__name__}: {e}"

# -- the graded run, timed end to end ----------------------------------------
t0 = time.perf_counter()
case = Case.from_yaml(case_path)
case.output.field = case.output.run_record = case.output.runtime_log = False
solver = Solver(case)
rec = solver.run()
wall_total = time.perf_counter() - t0

np.savez_compressed(out, vdf=np.asarray(solver.vdf), temp=np.asarray(rec.temp))

# -- the unit: one fine sweep *in this tree*, warm, median of five ----------
# The score is in the submission's own sweep units, so making the sweep
# faster buys nothing and only work outside the sweep is charged.  The
# verifier caps this from above by its own sweep time, so making it slower
# buys nothing either.  Part of the interface contract; absent => gate fails.
t_sweep = None
sweep_note = ""
try:
    solver.ctx.sweep(solver.mom, solver.vdf)                      # warm
    ts = []
    for _ in range(5):
        t0 = time.perf_counter()
        solver.ctx.sweep(solver.mom, solver.vdf)
        ts.append(time.perf_counter() - t0)
    t_sweep = float(np.median(ts))
except Exception as e:                                            # noqa: BLE001
    sweep_note = f"solver.ctx.sweep unavailable: {type(e).__name__}: {e}"

print(json.dumps({
    "iterations": int(rec.iterations),
    "converged": bool(rec.converged),
    "mass": float(rec.mass),
    "wall_total": wall_total,
    "t_sweep_agent_s": t_sweep,
    "sweep_note": sweep_note,
    "max_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
    "warm_up": warm_note or "ok",
}))
