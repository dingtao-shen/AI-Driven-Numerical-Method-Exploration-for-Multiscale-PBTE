"""Subprocess adapter for the frozen legacy runner and evaluator.

Only reviewed built-in code is executed. This is not an OS sandbox.
"""
from __future__ import annotations

import json
import runpy
import sys
from pathlib import Path

import numpy as np
import yaml


def main():
    mode, script, backend, case_path, dump, extra = sys.argv[1:]
    backend = Path(backend).resolve()
    sys.path.insert(0, str(backend))
    import pybte
    actual = Path(pybte.__file__).resolve()
    if actual.parent != backend / "pybte":
        raise RuntimeError(f"Wrong pybte import: {actual}")
    if mode == "candidate":
        spec = yaml.safe_load(Path(case_path).read_text())
        spec["scheme"] = {"method": "source" if extra == "cis" else "krylov", "accflag": 0}
        candidate_case = Path(dump).with_suffix(".candidate.yaml")
        candidate_case.write_text(yaml.safe_dump(spec, sort_keys=False))
        sys.argv = [script, str(backend), str(candidate_case), dump]
        result = runpy.run_path(script, run_name="__main__")
        rec = result["rec"]
        # The original runner dumps vdf before its sweep timing mutates the state.
        with np.load(dump, allow_pickle=False) as data:
            arrays = {k: data[k] for k in data.files}
        arrays.update(qx=rec.qx, qy=rec.qy, residual_history=rec.residual_history,
                      mass_history=rec.mass_history)
        np.savez_compressed(dump, **arrays)
        Path(dump).with_suffix(".details.json").write_text(json.dumps({
            "pybte_import": str(actual), "setup_wall_clock": rec.setup_wall_clock,
            "solver_wall_clock": rec.wall_clock, "factorisation_count": rec.factorisation_count,
            "layout": {"vdf": "direction, element, nodal DOF", "temp_qx_qy": "element integrals; no area normalization", "units": "legacy nondimensional Case flow normalization"},
        }, indent=2))
    elif mode == "evaluate":
        sys.argv = [script, str(backend), case_path, dump, extra]
        Path(dump).with_suffix(".verifier_import.json").write_text(json.dumps({"pybte_import": str(actual)}))
        runpy.run_path(script, run_name="__main__")
    else:
        raise ValueError(f"Unknown worker mode: {mode}")


if __name__ == "__main__":
    main()
