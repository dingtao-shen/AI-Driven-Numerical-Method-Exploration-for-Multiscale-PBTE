#!/usr/bin/env python3
"""Run k agent rollouts against a task and score each one.

This is the measurement Phase A exists for.  Everything before it -- the
ablation, the reference data, the oracle, the calibrated cells -- only
establishes that the task is well posed.  What the rollouts measure:

    gate pass rate   how often an agent delivers a solver that converges, to
                     the right fixed point, on every cell of every family;
    score            the speed-up over source iteration it achieved, which is
                     the benchmark's signal once gates are passed.  A task
                     that every agent passes is still informative if their
                     scores differ; one where every score is the same is not.

    0/k passing      read every transcript before believing it.  A task that
                     fails on import errors and unclear prompts is broken,
                     not hard.

Isolation
---------

``--isolation dir`` copies the environment to a scratch directory outside the
repository and runs the agent there.  Cheap, and enough to keep the agent from
stumbling over the answer -- but the full solver still exists elsewhere on the
filesystem, and an agent that goes looking with an absolute path could find
it.  Fine for a difficulty read, not for a published number.

``--isolation docker`` runs the agent inside a container that has the
environment and nothing else.  The solver, the oracle, the reference data and
the verifier are all outside it.  Use this for anything that goes in a paper.

In both modes the verifier runs afterwards, on the host, against the directory
the agent left behind.  It is never inside the sandbox: the transport-residual
check is only worth something if the submission never had access to the code
that computes it.

Usage::

    python tools/run_rollouts.py tasks/t01-square --k 5
    python tools/run_rollouts.py tasks/... --k 5 --isolation docker --model opus
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# The agent gets the task and a working copy of the solver.  It does not get
# the verifier, the reference data, the oracle, the ablation spec, or the
# unablated solver -- each of which contains the answer outright.
AGENT_SEES = ["environment", "task.md"]


def make_sandbox(task: Path, root: Path, trial: int) -> Path:
    box = root / f"trial_{trial:02d}"
    if box.exists():
        shutil.rmtree(box)
    box.mkdir(parents=True)
    for item in AGENT_SEES:
        src = task / item
        dst = box / item
        if src.is_dir():
            shutil.copytree(src, dst,
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        else:
            shutil.copy2(src, dst)
    return box


def run_agent_dir(box: Path, prompt: str, args) -> dict:
    cmd = [
        "claude", "-p", prompt,
        "--output-format", "stream-json", "--verbose",
        "--model", args.model,
        "--dangerously-skip-permissions",
    ]
    if args.max_turns:
        cmd += ["--max-turns", str(args.max_turns)]
    t0 = time.perf_counter()
    out = subprocess.run(cmd, cwd=box / "environment", capture_output=True,
                         text=True, timeout=args.timeout)
    return {"wall": time.perf_counter() - t0, "returncode": out.returncode,
            "stdout": out.stdout, "stderr": out.stderr[-4000:]}


def run_agent_docker(box: Path, prompt: str, args) -> dict:
    creds = Path.home() / ".claude"
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{box}:/task",
        "-v", f"{creds}:/root/.claude",
        "-w", "/task/environment",
        "-e", "OMP_NUM_THREADS=1", "-e", "MKL_NUM_THREADS=1",
        "-e", "OPENBLAS_NUM_THREADS=1", "-e", "NUMBA_NUM_THREADS=1",
        "-e", "PYTHONHASHSEED=0",
        args.image,
        "claude", "-p", prompt, "--output-format", "stream-json", "--verbose",
        "--model", args.model, "--dangerously-skip-permissions",
    ]
    if args.max_turns:
        cmd += ["--max-turns", str(args.max_turns)]
    t0 = time.perf_counter()
    out = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
    return {"wall": time.perf_counter() - t0, "returncode": out.returncode,
            "stdout": out.stdout, "stderr": out.stderr[-4000:]}


def parse_stream(stdout: str) -> dict:
    """The last ``result`` event of a stream-json transcript, or what there is.

    With ``--output-format json`` a run killed at the cap left nothing behind;
    the stream keeps every event up to the kill, and the final event -- when
    the run finishes -- carries turns, cost and duration.
    """
    meta = {"events": 0, "result_event": False}
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        meta["events"] += 1
        try:
            d = json.loads(line)
        except ValueError:
            continue
        if d.get("type") == "result":
            meta.update({k: d.get(k) for k in ("num_turns", "total_cost_usd", "duration_ms",
                                               "is_error", "subtype", "result",
                                               "api_error_status", "terminal_reason")})
            meta["result_event"] = True
    return meta


def voided_reason(meta: dict) -> str | None:
    """A trial the *harness* ended -- an API error such as a usage limit --
    measures nothing about the agent and must not count as a failure."""
    if meta.get("is_error") and meta.get("api_error_status"):
        return f"api error {meta['api_error_status']}: {str(meta.get('result'))[:80]}"
    return None


def archive_solution(task: Path, box: Path, out_dir: Path, trial: int) -> None:
    """Keep what the agent built, in the repository, at grading time.

    Sandboxes live under a scratch directory and do not survive a reboot; the
    scores and transcripts of a run are worthless for analysis without the
    code that produced them.  Two forms: a unified diff against the shipped
    environment (small, readable), and a tarball of the environment itself
    (complete, re-gradable under a later rule).
    """
    env = task / "environment"
    diff = subprocess.run(
        ["diff", "-ruN", "-x", "__pycache__", "-x", "*.pyc", "-x", ".pytest_cache",
         "-x", "out", str(env), str(box / "environment")],
        capture_output=True, text=True)
    (out_dir / f"trial_{trial:02d}_solution.diff").write_text(
        diff.stdout.replace(str(env), "environment").replace(str(box / "environment"), "solution"))
    import tarfile

    def _skip(ti):
        parts = Path(ti.name).parts
        return None if any(p in ("__pycache__", ".pytest_cache", "out") or p.endswith(".pyc")
                           for p in parts) else ti
    with tarfile.open(out_dir / f"trial_{trial:02d}_environment.tar.gz", "w:gz") as tf:
        tf.add(box / "environment", arcname="environment", filter=_skip)


def grade(task: Path, box: Path, out_json: Path) -> dict:
    r = subprocess.run(
        [sys.executable, str(task / "verifier" / "checks.py"),
         "--env", str(box / "environment"), "--json", str(out_json)],
        capture_output=True, text=True)
    if not out_json.exists():
        return {"gate": False, "score": None, "families": {},
                "verifier_error": (r.stderr or r.stdout)[-2000:]}
    return json.loads(out_json.read_text())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", type=Path)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--model", default="opus")
    ap.add_argument("--isolation", choices=["dir", "docker"], default="dir")
    ap.add_argument("--image", default="stiffkinetic-t01",
                    help="container image for --isolation docker")
    ap.add_argument("--timeout", type=float, default=7200.0,
                    help="per-rollout wall-clock cap, seconds")
    ap.add_argument("--max-turns", type=int, default=None)
    ap.add_argument("--sandbox-root", type=Path, default=None,
                    help="where to put the agents' working copies; must be "
                         "outside the repository")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--trials", type=int, nargs="*", default=None,
                    help="trial numbers to run (default 1..k); use to replace "
                         "voided trials -- earlier rows in summary.json are kept")
    args = ap.parse_args(argv)

    task = args.task.resolve()
    prompt = (task / "task.md").read_text()
    out_dir = args.out or (ROOT / "experiments" / "results" / task.name)
    out_dir.mkdir(parents=True, exist_ok=True)
    sandbox_root = args.sandbox_root or Path(tempfile.mkdtemp(prefix="rollout-"))
    if ROOT in sandbox_root.resolve().parents or sandbox_root.resolve() == ROOT:
        print("refusing to put sandboxes inside the repository -- the agent "
              "would be one `cd ..` away from the solver", file=sys.stderr)
        return 2

    print(f"task       {task.name}")
    print(f"model      {args.model}   isolation {args.isolation}   k={args.k}")
    print(f"sandboxes  {sandbox_root}")
    print(f"results    {out_dir}\n")

    runner = run_agent_docker if args.isolation == "docker" else run_agent_dir
    trials = args.trials or list(range(1, args.k + 1))
    rows = []
    for i in trials:
        box = make_sandbox(task, sandbox_root, i)
        print(f"trial {i}: agent running ...", flush=True)
        try:
            agent = runner(box, prompt, args)
        except subprocess.TimeoutExpired as e:
            # The sandbox is graded as the agent left it -- a solution that
            # was in place before the cap still counts -- but the transcript
            # is whatever had been flushed.  With --output-format json that
            # is usually nothing; stream-json would keep the events.
            def _txt(b):
                return b.decode(errors="replace") if isinstance(b, bytes) else (b or "")
            agent = {"wall": args.timeout, "returncode": -1, "timed_out": True,
                     "stdout": _txt(e.stdout), "stderr": "timeout\n" + _txt(e.stderr)[-4000:]}
        (out_dir / f"trial_{i:02d}_transcript.jsonl").write_text(agent["stdout"] or "")
        meta = parse_stream(agent["stdout"] or "")
        meta["timed_out"] = bool(agent.get("timed_out"))
        (out_dir / f"trial_{i:02d}_meta.json").write_text(json.dumps(meta, indent=2))
        if agent["stderr"]:
            (out_dir / f"trial_{i:02d}_stderr.txt").write_text(agent["stderr"])

        archive_solution(task, box, out_dir, i)
        res = grade(task, box, out_dir / f"trial_{i:02d}_score.json")
        cells = [c for f in res.get("families", {}).values() for c in f["cells"]]
        n = sum(1 for c in cells if c.get("passed"))
        fam = {k: (f["gate"], f["score"]) for k, f in res.get("families", {}).items()}
        voided = voided_reason(meta)
        rows.append({"trial": i, "gate": bool(res.get("gate")), "score": res.get("score"),
                     "families": fam, "cells_passed": n, "cells": len(cells),
                     "agent_wall": agent["wall"], "timed_out": bool(agent.get("timed_out")),
                     "turns": meta.get("num_turns"), "cost_usd": meta.get("total_cost_usd"),
                     "events": meta.get("events"),
                     "agent_returncode": agent["returncode"],
                     "voided": voided,
                     "sandbox": str(box)})
        sc = f"{res['score']:.1f}x" if res.get("score") else "-"
        print(f"trial {i}: gate {'PASS' if res.get('gate') else 'FAIL'}   score {sc}   "
              f"{n}/{len(cells)} cells   agent {agent['wall'] / 60:.0f} min"
              f"{'   (timed out; graded as left)' if agent.get('timed_out') else ''}"
              f"{f'   VOIDED -- {voided}' if voided else ''}\n", flush=True)

    # Rows from an earlier invocation (a run resumed, or voided trials being
    # replaced with --trials) are kept; a re-run trial number replaces its row.
    summary_path = out_dir / "summary.json"
    if summary_path.exists():
        old = json.loads(summary_path.read_text()).get("trials", [])
        done = {r["trial"] for r in rows}
        rows = sorted(rows + [r for r in old if r["trial"] not in done],
                      key=lambda r: r["trial"])
    counted = [r for r in rows if not r.get("voided")]
    solved = sum(1 for r in counted if r["gate"])
    summary = {"task": task.name, "model": args.model, "k": args.k,
               "isolation": args.isolation, "solved": solved,
               "voided": len(rows) - len(counted), "trials": rows}
    summary_path.write_text(json.dumps(summary, indent=2))

    print(f"{'trial':>6} {'gate':>5} {'score':>8} {'cells':>6} {'minutes':>8}")
    for r in rows:
        sc = f"{r['score']:.1f}x" if r["score"] else "-"
        print(f"{r['trial']:6d} {str(r['gate']):>5} {sc:>8} "
              f"{r['cells_passed']:3d}/{r['cells']:<2d} {r['agent_wall'] / 60:8.0f}"
              f"{'   voided' if r.get('voided') else ''}")
    print(f"\n{solved}/{len(counted)} passed every gate"
          + (f"   ({len(rows) - len(counted)} voided, not counted)" if len(rows) != len(counted) else ""))
    if solved == 0:
        print("0/k -- read every transcript before concluding the task is hard; "
              "failures on imports, paths or an unclear prompt mean it is broken")
    scores = [r["score"] for r in rows if r["score"]]
    if len(scores) >= 2:
        print(f"score range {min(scores):.1f}x .. {max(scores):.1f}x -- "
              "spread across trials/models is the benchmark's signal now, not pass rate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
