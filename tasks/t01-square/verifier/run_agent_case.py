"""Run one case in the *agent's* tree and dump what the verifier needs.

Touches only the public interface the task promises -- ``Case.from_yaml``,
``Solver(case).run()``, ``solver.vdf`` -- so a legitimate rewrite of the
internals cannot break it.

Timing is the whole of ``Case.from_yaml + Solver(case) + run()`` for the
graded case, after a warm-up solve on a *different mesh*: the warm-up gets
kernels compiled and caches loaded so the measurement is of the method and
not of the JIT, while a different mesh means nothing about the graded case
can have been precomputed or memoised.  Setup is inside the timing on
purpose -- a method that moves all its work into "setup" still pays for it.
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

# -- warm-up on another mesh, all walls thermalising, three iterations ----
warm_note = ""
try:
    w = copy.deepcopy(spec)
    w["mesh"]["file"] = str(Path(spec["mesh"]["file"]).parent / "C1_Nx11_Ny11.msh")
    w["flow"]["tau_r"] = 1.0
    w["iteration"]["tmax"] = 3
    for b in w["boundaries"]:
        b["type"] = "thermalising"
        b.pop("xoff", None)
        b.pop("yoff", None)
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
print(json.dumps({
    "iterations": int(rec.iterations),
    "converged": bool(rec.converged),
    "mass": float(rec.mass),
    "wall_total": wall_total,
    "max_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
    "warm_up": warm_note or "ok",
}))
