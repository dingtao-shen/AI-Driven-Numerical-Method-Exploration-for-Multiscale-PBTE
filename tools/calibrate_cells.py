#!/usr/bin/env python3
"""Calibrate a task's cells: the unaccelerated cost, and the fixed-point tolerance.

For every case the task ships this measures, on the shipped discretisation:

* the unaccelerated solver's iteration count to ``tol`` (or that it did not
  converge within ``--cis-cap``, in which case the count is a censored lower
  bound) and its wall-clock per iteration on this machine;
* the oracle's iteration count, and -- as a diagnostic only -- the spread
  between the oracle stopped at the shipped ``tol`` and two decades tighter.
  The Krylov oracle converges to round-off either way, so that spread is
  ~0 and cannot set a threshold.  The field gate is an explicit, uniform
  ``--field-gap-max`` (default 1e-5): an order above the reference's own
  accuracy, two orders below any method on a displaced fixed point.

A cell is *scored* when the unaccelerated solver needs at least
``--slow-from`` iterations there; elsewhere it is already cheap and the cell
is a gate only.  No budget is derived from anything: the score is a speed-up,
which needs the unaccelerated count and nothing else.

Usage::

    python tools/calibrate_cells.py tasks/t01-square
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent


def round_up_1sf(x: float) -> float:
    if x <= 0:
        return 0.0
    e = math.floor(math.log10(x))
    return math.ceil(x / 10 ** e) * 10 ** e


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", type=Path)
    ap.add_argument("--cis-cap", type=int, default=200_000)
    ap.add_argument("--slow-from", type=int, default=1000)
    ap.add_argument("--transport-residual-max", type=float, default=1e-7)
    ap.add_argument("--field-gap-max", type=float, default=1e-5)
    args = ap.parse_args(argv)

    task = args.task.resolve()
    spec = yaml.safe_load((task / "ablation.yaml").read_text())
    src = (task / spec["source"]).resolve()
    sys.path.insert(0, str(src))
    from pybte import Case, Solver                                # noqa: E402

    cases = {k[len("cases/"):-len(".yaml")]: v
             for k, v in spec["new_files"].items() if k.startswith("cases/")}
    fam_desc = {}
    for name, text in cases.items():
        fam_desc.setdefault(name.split("_")[0], text.splitlines()[0].lstrip("# ").split("  Kn_R")[0])

    families: dict = {}
    print(f"{'cell':>12} | {'CIS it':>8} {'s/it':>7} | {'oracle':>7} {'spread':>9} "
          f"{'gate':>7} | scored")
    print("-" * 72)
    with tempfile.TemporaryDirectory() as td:
        for name, text in cases.items():
            fam = name.split("_")[0]
            spec_y = yaml.safe_load(text)
            spec_y["mesh"]["file"] = str(src / "meshes" / Path(spec_y["mesh"]["file"]).name)
            cf = Path(td) / f"{name}.yaml"
            cf.write_text(yaml.safe_dump(spec_y, sort_keys=False))

            def run(method, tol, tmax):
                c = Case.from_yaml(cf)
                c.scheme.accflag = 0
                c.scheme.method = method
                c.iteration.tol = tol
                c.iteration.tmax = tmax
                c.output.field = c.output.run_record = c.output.runtime_log = False
                t0 = time.perf_counter()
                s = Solver(c)
                r = s.run()
                return r, time.perf_counter() - t0

            rc, wc = run("source", spec_y["iteration"]["tol"], args.cis_cap)
            ro, _ = run("krylov", spec_y["iteration"]["tol"], args.cis_cap)
            rt, _ = run("krylov", spec_y["iteration"]["tol"] * 1e-2, args.cis_cap)
            spread = float(np.abs(ro.temp - rt.temp).max() / np.abs(rt.temp).max())
            gate = args.field_gap_max
            cis_it = int(rc.iterations) if rc.converged else None
            scored = cis_it is None or cis_it >= args.slow_from
            cell = {
                "kn_r": float(spec_y["flow"]["tau_r"]), "kn_n": float(spec_y["flow"]["tau_n"]),
                "cis_iterations": cis_it, "cis_cap": args.cis_cap,
                "cis_seconds_per_iteration": wc / max(1, rc.iterations),
                "oracle_iterations": int(ro.iterations),
                "oracle_transport_residual": float(ro.diagnostics["transport_residual"]),
                "measured_spread": spread, "field_gap_max": gate, "scored": scored,
            }
            families.setdefault(fam, {"description": fam_desc[fam], "cells": {}})
            families[fam]["cells"][name] = cell
            print(f"{name:>12} | {(str(cis_it) if cis_it else f'>{args.cis_cap}'):>8} "
                  f"{cell['cis_seconds_per_iteration']:7.4f} | {ro.iterations:7d} "
                  f"{spread:9.2e} {gate:7.0e} | {scored}", flush=True)

    out = task / "verifier" / "cells.json"
    out.write_text(json.dumps({
        "_comment": "Generated by tools/calibrate_cells.py; see its docstring for the rules.",
        "transport_residual_max": args.transport_residual_max,
        "field_gap_max": args.field_gap_max,
        "families": families}, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
