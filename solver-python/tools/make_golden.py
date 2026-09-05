#!/usr/bin/env python3
"""Regenerate everything under ``fortran-reference/golden/``.

The golden data is ~600 MB of raw stage dumps and residual histories, so it is
not committed.  This script rebuilds it from the reference sources.  Run it
once after cloning if you want the ``fortran`` and ``stage`` tests to do
anything other than skip.

    solver-python/tools/build_fortran.sh --variant A
    solver-python/tools/build_fortran.sh --variant B
    python solver-python/tools/make_golden.py [--quick] [--only dumps|runs]

``--quick`` skips the two multi-hour CIS runs at ``tau_R <= 1e-2``.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REPO = ROOT.parent
GOLDEN = REPO / "fortran-reference" / "golden"
RUNNER = HERE / "run_fortran.py"

#: stage dumps: (name, extra args).  Cheap -- 10 iterations each.
DUMPS = [
    ("dump_cis_shipped",   ["--accflag", "0"]),
    ("dump_cis_deg1",      ["--accflag", "0", "--deg", "1"]),
    ("dump_cis_deg2",      ["--accflag", "0", "--deg", "2"]),
    ("dump_cis_np10na20",  ["--accflag", "0", "--npole", "10", "--nazim", "20"]),
    ("dump_cis_tauN1e-2",  ["--accflag", "0", "--tau-n", "1e-2"]),
    ("dump_gsis_shipped",  ["--accflag", "1"]),
    ("dump_gsis_deg1",     ["--accflag", "1", "--deg", "1"]),
    ("dump_gsisB_shipped", ["--accflag", "1", "--variant", "B"]),
]

#: full runs, for the integration tests.  (name, args, minutes, quick?)
RUNS = [
    ("cis_tauR1.0_deg3",   ["--accflag", "0", "--tau-r", "1.0"], 0.1, True),
    ("cis_tauR1e-1_deg3",  ["--accflag", "0", "--tau-r", "1e-1"], 1, True),
    ("gsis_tauR1e-3_deg3", ["--accflag", "1", "--tau-r", "1e-3"], 0.2, True),
    ("fx_cis_1e-1",        ["--accflag", "0", "--tau-r", "1e-1",
                            "--dump", "--dump-iters", "1", "--dump-no-vdf"], 1, True),
    ("fx_gsis_1e-1",       ["--accflag", "1", "--tau-r", "1e-1",
                            "--dump", "--dump-iters", "1", "--dump-no-vdf"], 0.2, True),
    ("cis_tauR1e-2_deg3",  ["--accflag", "0", "--tau-r", "1e-2"], 58, False),
    ("cis_tauR1e-3_deg3_trunc20000",
                           ["--accflag", "0", "--tau-r", "1e-3",
                            "--tmax", "20000"], 68, False),
]


def run(name, extra, dump=False):
    out = GOLDEN / name
    cmd = [sys.executable, str(RUNNER), "--out", str(out), "--tmax", "8000000"]
    if dump:
        cmd += ["--dump", "--dump-iters", "1,2,3,10", "--tmax", "10"]
    cmd += extra
    t0 = time.time()
    rc = subprocess.run(cmd).returncode
    print(f"  {'ok ' if rc == 0 else 'FAIL'} {name}  ({time.time() - t0:.0f}s)")
    return rc


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--quick", action="store_true",
                   help="skip the two multi-hour CIS runs")
    p.add_argument("--only", choices=["dumps", "runs"], default=None)
    args = p.parse_args(argv)

    for flavour in ("varA", "varB"):
        exe = REPO / "fortran-build" / flavour / "DGACC"
        if not exe.exists():
            print(f"missing {exe}\n  run: solver-python/tools/build_fortran.sh "
                  f"--variant {flavour[-1]}", file=sys.stderr)
            return 2

    GOLDEN.mkdir(parents=True, exist_ok=True)
    bad = 0

    if args.only != "runs":
        print("== stage dumps ==")
        for name, extra in DUMPS:
            variant = ["--variant", "B"] if "gsisB" in name else []
            bad += bool(run(name, extra + variant, dump=True))

    if args.only != "dumps":
        print("== full runs ==")
        skipped = []
        for name, extra, minutes, quick_ok in RUNS:
            if args.quick and not quick_ok:
                skipped.append((name, minutes))
                continue
            print(f"  ... {name} (about {minutes:g} min)")
            bad += bool(run(name, extra))
        for name, minutes in skipped:
            print(f"  -- {name} skipped (--quick; about {minutes:g} min)")

    print(f"\n{'all good' if not bad else f'{bad} failures'} -> {GOLDEN}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
