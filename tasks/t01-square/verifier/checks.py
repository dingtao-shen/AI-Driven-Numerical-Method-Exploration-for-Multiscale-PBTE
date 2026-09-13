#!/usr/bin/env python3
"""Score a submission: gates per cell, a score per family, one overall.

The design separates what is non-negotiable from what is a matter of degree.

**Gates** (pass/fail, every cell of every family):

1. *converged* -- the run reports convergence inside the iteration cap.
2. *transport residual* -- ``||A f - b|| / ||b||`` of the discrete kinetic
   system, computed from the submitted distribution **in a pristine tree the
   submission never touched**.  A state at round-off transport residual is
   the discrete fixed point; loosening a tolerance or damping the physics
   cannot move this number.
3. *field* -- the converged temperature matches the certified reference,
   both as reported and as recomputed from the submitted distribution, to
   ``field_gap_max`` (1e-5, L-inf relative): an order above the reference's
   own accuracy, two orders below any method on a displaced fixed point.
   Stopping when successive iterates differ by 1e-8 does not meet it at
   small Knudsen number -- that is the pseudo-convergence the benchmark is
   about, and it is caught here as well as by gate 2.

A family scores only if every one of its cells passes every gate; the
benchmark scores only if every family does.

4. *interface* -- ``solver.ctx.sweep(solver.mom, solver.vdf)`` is callable
   and performs one fine transport sweep; the verifier times it to set the
   unit below.

Every cell is run in a *copy* of the submission with no compiled-kernel
cache, after an untimed warm-up of the same kind of case on a coarser mesh
(``run_agent_case.py``), so the timing is of the method and is reproducible
from the archived tree.

**Score**: on the cells where the unaccelerated solver is slow, the speed-up
over it in *sweep-equivalents*, combined as a geometric mean.  The
submission's cost is its ``Case + Solver + run`` wall-clock (setup included)
divided by **its own** fine-sweep time, both timed by the verifier in the
same single-threaded session -- with the unit capped from above by the
verifier's own sweep time.  Making the sweep faster shrinks the unit with
the wall-clock and gains nothing; making it slower is charged in the
verifier's unit; only work outside the sweep -- a low-order solve, a
factorisation -- moves the number, and it is paid for at its actual cost in
sweep units.  The unaccelerated side is its setup plus its calibrated
iteration count times a freshly timed seconds-per-iteration, in the
verifier's sweep units.  The self-reported sweep count is shown alongside
but does not enter the score.  Where the unaccelerated solver never
converged within its calibration cap the speed-up is a lower bound.

Usage::

    python verifier/checks.py [--env environment] [--pristine ...] [--json out]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
TASK = HERE.parent


# Grading is single-threaded, on both sides of the ratio.  The score is a
# ratio of wall-clocks; letting either side use however many cores it finds
# would make it a statement about the machine, and a submission that
# parallelises the sweep would be scoring core count rather than iteration
# design.  (The Dockerfile pins the same variables for container runs.)
PINNED = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
              OPENBLAS_NUM_THREADS="1", NUMBA_NUM_THREADS="1",
              PYTHONHASHSEED="0")


def run_stage(script, *args):
    out = subprocess.run([sys.executable, str(HERE / script), *map(str, args)],
                         capture_output=True, text=True, env=PINNED)
    if out.returncode != 0:
        return None, (out.stderr or out.stdout).strip()[-2000:]
    try:
        return json.loads(out.stdout.strip().splitlines()[-1]), None
    except (ValueError, IndexError):
        return None, f"stage {script} produced no result:\n{out.stdout[-2000:]}"


def stage_case(name, pristine, tmp):
    """The verifier grades on *its own* copy of the case and mesh."""
    spec = yaml.safe_load((HERE / "cases" / f"{name}.yaml").read_text())
    spec["mesh"]["file"] = str(pristine / "meshes" / Path(spec["mesh"]["file"]).name)
    out = tmp / f"{name}.yaml"
    out.write_text(yaml.safe_dump(spec, sort_keys=False))
    return out


def score_cell(name, cell, env, pristine, spec, tmp):
    res = {"cell": name, "gates": {}}
    case_path = stage_case(name, pristine, tmp)
    ref = TASK / "reference" / f"{name}.npz"
    dump = tmp / f"{name}.npz"

    run, err = run_stage("run_agent_case.py", env, case_path, dump)
    if err:
        res.update(error=err, passed=False)
        return res
    if run.get("warm_up") != "ok":
        # A cold run is charged for the JIT; that is a verifier fault, not a
        # property of the submission, and the cell must not be scored on it.
        res.update(run, error=f"verifier warm-up did not run: {run.get('warm_up')}",
                   passed=False)
        return res
    res.update(run)
    ev, err = run_stage("evaluate.py", pristine, case_path, dump, ref)
    if err or "error" in (ev or {}):
        res.update(error=err or ev["error"], passed=False)
        return res
    res.update(ev)

    g = res["gates"]
    g["converged"] = bool(run["converged"])
    g["transport_residual"] = ev["transport_residual"] < spec["transport_residual_max"]
    gap = max(ev["field_gap"], ev["field_gap_from_vdf"])
    res["field_gap_worst"] = gap
    res["field_gap_max"] = cell["field_gap_max"]
    g["field"] = gap < cell["field_gap_max"]
    g["interface"] = run.get("t_sweep_agent_s") is not None      # ctx.sweep callable
    res["passed"] = all(g.values())

    # Cost in the submission's own sweep units, capped by the verifier's:
    #   E_agent = wall / min(own sweep, pristine sweep)
    # A faster sweep shrinks the unit with the wall-clock and gains nothing;
    # a slower one is charged in the pristine unit.  Only work outside the
    # sweep changes E_agent.  The unaccelerated side is in pristine units.
    t_p = ev["t_sweep_s"]
    t_a = run.get("t_sweep_agent_s") or t_p
    unit = min(t_a, t_p)
    res["t_sweep_agent_s"] = run.get("t_sweep_agent_s")
    res["t_sweep_pristine_s"] = t_p
    res["sweep_equivalents"] = run["wall_total"] / unit
    n_cis = cell.get("cis_iterations")
    censored = n_cis is None
    n_cis = cell["cis_cap"] if censored else n_cis
    res["cis_sweep_equivalents"] = (ev["cis_setup_s"] + n_cis * ev["cis_s_per_iter"]) / t_p
    res["scored"] = bool(cell.get("scored", False))
    if res["scored"]:
        res["speedup"] = res["cis_sweep_equivalents"] / max(res["sweep_equivalents"], 1e-9)
        res["speedup_censored"] = censored
    return res


def geo_mean(xs):
    xs = [x for x in xs if x and x > 0]
    return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--env", type=Path, default=TASK / "environment")
    ap.add_argument("--pristine", type=Path, default=TASK / "verifier" / "pristine")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--families", nargs="*", default=None,
                    help="restrict to these families (default: all)")
    args = ap.parse_args(argv)

    spec = json.loads((HERE / "cells.json").read_text())
    env = args.env.resolve()
    pristine = args.pristine.resolve()
    if not (pristine / "pybte").is_dir():
        print(f"verifier: no pristine solver at {pristine}", file=sys.stderr)
        return 2

    families = {}
    with tempfile.TemporaryDirectory() as td:
        # Grade a copy of the submission with no compiled-kernel cache in it:
        # what a submission happened to run before grading must not change
        # its timing, and a re-score of the archived tree must reproduce it.
        graded = Path(td) / "environment"
        shutil.copytree(env, graded, symlinks=True,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".pytest_cache", "out"))
        env = graded
        for fam, fspec in spec["families"].items():
            if args.families and fam not in args.families:
                continue
            cells = [score_cell(n, c, env, pristine, spec, Path(td))
                     for n, c in fspec["cells"].items()]
            gate = all(c.get("passed") for c in cells)
            scored = [c for c in cells if c.get("scored") and "speedup" in c]
            score = geo_mean([c["speedup"] for c in scored]) if gate else None
            families[fam] = {
                "description": fspec.get("description", ""),
                "gate": gate, "score": score,
                "score_censored": any(c.get("speedup_censored") for c in scored),
                "cells": cells,
                "max_rss_mb": max((c.get("max_rss_mb") or 0) for c in cells) or None,
            }

    overall_gate = all(f["gate"] for f in families.values())
    overall = geo_mean([f["score"] for f in families.values()]) if overall_gate else None

    hdr = (f"{'cell':>12} {'sweeps':>7} {'sweep-eq':>9} {'conv':>5} {'transport':>10} "
           f"{'field gap':>10} {'gates':>6} {'speed-up':>10}")
    for fam, f in families.items():
        print(f"\n[{fam}] {f['description']}")
        print(hdr)
        print("-" * len(hdr))
        for c in f["cells"]:
            if "error" in c:
                print(f"{c['cell']:>12}  ERROR: {c['error'].splitlines()[-1][:60]}")
                continue
            marks = "".join("." if v else "X" for v in c["gates"].values())
            sp = ""
            if c.get("scored"):
                sp = f"{'>=' if c.get('speedup_censored') else ''}{c['speedup']:.1f}x"
            print(f"{c['cell']:>12} {c['iterations']:7d} {c['sweep_equivalents']:9.0f} "
                  f"{str(c['gates']['converged']):>5} "
                  f"{c['transport_residual']:10.2e} {c['field_gap_worst']:10.2e} "
                  f"{marks:>6} {sp:>10}")
        line = (f"  gate {'PASS' if f['gate'] else 'FAIL'}")
        if f["score"] is not None:
            line += (f"   score {'>=' if f['score_censored'] else ''}{f['score']:.1f}x")
        if f["max_rss_mb"]:
            line += f"   peak RSS {f['max_rss_mb']:.0f} MB"
        print(line)

    print(f"\noverall gate {'PASS' if overall_gate else 'FAIL'}"
          + (f"   score {overall:.1f}x" if overall is not None else "   (no score)"))
    if args.json:
        args.json.write_text(json.dumps(
            {"gate": overall_gate, "score": overall, "families": families}, indent=2))
    return 0 if overall_gate else 1


if __name__ == "__main__":
    raise SystemExit(main())
