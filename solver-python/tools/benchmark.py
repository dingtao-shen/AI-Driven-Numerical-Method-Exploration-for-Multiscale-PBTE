#!/usr/bin/env python3
"""Measure the performance targets.

    python tools/benchmark.py [--json bench.json]

Reports setup cost, per-iteration cost and peak RSS for the shipped case, and
checks each against its target.  Run it on an otherwise idle
machine: the numbers are wall-clock.
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from pybte import Case, Solver                                    # noqa: E402

CASE = ROOT / "cases" / "cavity_tauR1e-3_cis.yaml"

#: (label, target, unit, comparison)
TARGETS = [
    ("setup", 30.0, "s"),
    ("gsis_setup", 30.0, "s"),
    ("cis_per_iter_ms", 50.0, "ms"),
    ("gsis_per_iter_ms", 150.0, "ms"),
    ("gsis_full_tauR1e-3", 60.0, "s"),
    ("cis_full_tauR1e-1", 300.0, "s"),
    ("peak_rss_mb", 2048.0, "MB"),
]


def _case(accflag=0, tau_r=1e-3, tmax=8_000_000):
    c = Case.from_yaml(CASE)
    c.scheme.accflag = accflag
    c.flow.tau_r = tau_r
    c.iteration.tmax = tmax
    c.output.field = c.output.run_record = c.output.runtime_log = False
    return c


def peak_rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def time_per_iteration(solver, n=40, warmup=3):
    for _ in range(warmup):
        solver.step()
    t0 = time.perf_counter()
    for _ in range(n):
        solver.step()
    return (time.perf_counter() - t0) / n * 1e3


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--json", default=None)
    args = p.parse_args(argv)

    out = {}

    # everything below the first sweep is setup: mesh, basis, integration,
    # ordering and the operator precompute
    t0 = time.perf_counter()
    s = Solver(_case(0))
    out["setup"] = time.perf_counter() - t0
    out["sweep_mode"] = s.ctx.mode
    out["operator_mb"] = s.ctx.lu.nbytes / 1e6 if s.ctx.lu is not None else 0.0
    out["n_tris"] = s.n_tris
    out["ndir"] = s.ndir
    out["cis_per_iter_ms"] = time_per_iteration(s)

    t0 = time.perf_counter()
    sg = Solver(_case(1))
    out["gsis_setup"] = time.perf_counter() - t0
    out["global_matrix_n"] = int(sg.acc.K.shape[0])
    out["global_matrix_nnz"] = int(sg.acc.K.nnz)
    out["gsis_per_iter_ms"] = time_per_iteration(sg)

    t0 = time.perf_counter()
    r = Solver(_case(1, 1e-3)).run()
    out["gsis_full_tauR1e-3"] = time.perf_counter() - t0
    out["gsis_full_iters"] = r.iterations

    t0 = time.perf_counter()
    r = Solver(_case(0, 1e-1)).run()
    out["cis_full_tauR1e-1"] = time.perf_counter() - t0
    out["cis_full_iters"] = r.iterations

    out["peak_rss_mb"] = peak_rss_mb()

    print(f"shipped case: N_TRIS={out['n_tris']}, {out['ndir']} directions, "
          f"DEG=3, sweep mode {out['sweep_mode']} ({out['operator_mb']:.0f} MB)")
    print(f"GSIS global matrix: {out['global_matrix_n']} x "
          f"{out['global_matrix_n']}, {out['global_matrix_nnz']} nonzeros\n")
    print(f"{'quantity':<26} {'measured':>12} {'target':>10}   verdict")
    print("-" * 62)
    ok = True
    for key, target, unit in TARGETS:
        v = out[key]
        good = v < target
        ok &= good
        print(f"{key:<26} {v:>9.2f} {unit:<2} {target:>7.0f} {unit:<2}  "
              f"{'PASS' if good else 'FAIL'}  ({target / v:.0f}x margin)"
              if good else
              f"{key:<26} {v:>9.2f} {unit:<2} {target:>7.0f} {unit:<2}  FAIL")
    print(f"\nGSIS at tau_R=1e-3 converged in {out['gsis_full_iters']} iterations, "
          f"CIS at tau_R=1e-1 in {out['cis_full_iters']}.")

    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=2))
        print(f"wrote {args.json}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
