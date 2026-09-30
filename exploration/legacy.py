"""Compatibility glue; legacy formulas and gates remain in frozen checks.py."""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from . import BASELINE, ROOT
from .contracts import read_yaml
from .records import clean_copy, sha256, verify_files, write_json

THREADS = {"OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
           "NUMBA_NUM_THREADS": "1", "PYTHONHASHSEED": "0"}
TASK = ROOT / "tasks/t01-square"
WORKER = ROOT / "exploration/worker.py"


def load_profile(path: Path) -> dict:
    profile = read_yaml(path)
    required = {"schema_version": 1, "profile_id": "legacy_t01_smoke",
                "purpose": "compatibility_smoke", "protocol_id": "legacy_t01_v0_3_2",
                "counts_as_new_research_evidence": False, "live_agent_enabled": False,
                "source_revision": "dffff2319b290860b5c4a87509428595bce4d2b0"}
    for key, value in required.items():
        if profile.get(key) != value or type(profile.get(key)) is not type(value):
            raise ValueError(f"Unknown or invalid legacy profile field: {key}")
    cases = profile.get("cases")
    if not isinstance(cases, list) or not cases or len(set(cases)) != len(cases):
        raise ValueError("profile cases must be a nonempty unique list")
    if any(c not in ("F1_1_1", "F2_1_1", "F3_1_1") for c in cases):
        raise ValueError("This engineering profile only admits existing F1/F2/F3 *_1_1 cases")
    manifest = json.loads(BASELINE.read_text())
    if profile.get("verifier_revision") != manifest["provenance"]["verifier_revision"]:
        raise ValueError("Profile verifier revision mismatch")
    if profile.get("baseline_manifest_sha256") != sha256(BASELINE):
        raise ValueError("Profile baseline manifest SHA256 mismatch")
    files = {p: h for p, h in manifest["files"].items() if p.startswith(
        ("tasks/t01-square/", "solver-python/pybte/"))}
    issues = verify_files(files)
    if issues:
        raise ValueError("; ".join(issues))
    return profile


def prepare_trusted(out: Path) -> Path:
    """Freeze trusted fixture separately; Python caches live only in the run."""
    target = out / "trusted/t01-square"
    target.mkdir(parents=True, exist_ok=False)
    baseline = json.loads(BASELINE.read_text())
    files = {p: h for p, h in baseline["files"].items() if p.startswith("tasks/t01-square/")}
    issues = verify_files(files)
    if issues:
        raise ValueError("; ".join(issues))
    # Ignore any personal/untracked file that happens to be in the fixture tree.
    for name in files:
        path = target / Path(name).relative_to("tasks/t01-square")
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, path)
        path.chmod(0o444)
    return target


def validate_solution(path: Path) -> None:
    import numpy as np
    with np.load(path, allow_pickle=False) as data:
        for key in ("vdf", "temp"):
            value = data[key]
            if not np.issubdtype(value.dtype, np.number) or not np.isfinite(value).all():
                raise ValueError(f"non-finite or nonnumeric {key}")
        # Exact shape checks are performed by the pristine evaluator.
        if data["vdf"].ndim != 3 or data["temp"].ndim != 1:
            raise ValueError("Invalid solution dimensions")


def run_stage(command: list[str], out: Path, name: str, timeout: float, env: dict) -> tuple[dict | None, str | None, dict]:
    from .records import now
    meta = {"started_at": now(), "command": command, "timeout_s": timeout}
    try:
        result = subprocess.run(command, cwd=out, env=env, capture_output=True,
                                text=True, timeout=timeout)
        stdout, stderr = result.stdout, result.stderr
        meta.update(returncode=result.returncode, termination="completed" if result.returncode == 0 else "execution_error")
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or b""
        stderr = exc.stderr or b""
        stdout = stdout.decode(errors="replace") if isinstance(stdout, bytes) else stdout
        stderr = stderr.decode(errors="replace") if isinstance(stderr, bytes) else stderr
        meta.update(returncode=None, termination="timeout")
    meta["finished_at"] = now()
    (out / f"{name}.stdout.log").write_text(stdout)
    (out / f"{name}.stderr.log").write_text(stderr)
    if meta["termination"] != "completed":
        write_json(out / f"{name}.process.json", meta)
        return None, f"{name}: {meta['termination']}; {(stderr or stdout)[-1000:]}", meta
    try:
        data = json.loads(stdout.strip().splitlines()[-1], parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
        if not isinstance(data, dict):
            raise ValueError("Expected result mapping")
    except (ValueError, IndexError):
        meta["termination"] = "invalid_output"
        write_json(out / f"{name}.process.json", meta)
        return None, f"{name}: invalid or missing JSON result", meta
    write_json(out / f"{name}.process.json", meta)
    return data, None, meta


def evaluate_cell(trusted: Path, candidate: Path, name: str, method: str,
                  out: Path, timeout: float = 180) -> tuple[dict, list[dict]]:
    """Call the original score_cell; replace only subprocess/file boundaries."""
    spec = importlib.util.spec_from_file_location("legacy_checks", trusted / "verifier/checks.py")
    checks = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checks)
    config = json.loads((trusted / "verifier/cells.json").read_text())
    cell = config["families"][name.split("_")[0]]["cells"][name]
    pristine = trusted / "verifier/pristine"
    work = out / "candidate_work"
    clean_copy(candidate, work)  # per-cell cache-free tree, same as the old checker
    processes = []
    env = {k: os.environ[k] for k in ("PATH", "LANG", "LC_ALL", "SYSTEMROOT") if k in os.environ}
    env.update(THREADS, PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1",
               NUMBA_CACHE_DIR=str(out / "numba_cache"), TMPDIR=str(out))

    def stage(script, *args):
        if script == "run_agent_case.py":
            backend, case, dump = args
            command = [sys.executable, str(WORKER), "candidate", str(checks.HERE / script),
                       str(backend), str(case), str(dump), method]
            tag = "candidate"
        else:
            backend, case, dump, reference = args
            try:
                validate_solution(Path(dump))
            except (ValueError, KeyError, OSError) as exc:
                processes.append({"termination": "invalid_output", "error": str(exc)})
                return None, str(exc)
            command = [sys.executable, str(WORKER), "evaluate", str(checks.HERE / script),
                       str(backend), str(case), str(dump), str(reference)]
            tag = "verifier"
            # A distinct cache isolates trusted JIT from candidate modules.
        local_env = dict(env, NUMBA_CACHE_DIR=str(out / f"{tag}_numba_cache"))
        data, error, meta = run_stage(command, out, tag, timeout, local_env)
        processes.append(meta)
        return data, error

    checks.run_stage = stage
    result = checks.score_cell(name, cell, work, pristine, config, out)
    # Gate formulas above are unchanged. No new-protocol score is derived here.
    write_json(out / "raw_legacy.json", result)
    return result, processes
