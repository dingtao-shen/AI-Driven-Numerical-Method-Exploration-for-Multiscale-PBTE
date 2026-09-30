"""Run, independently evaluate and archive a controlled compatibility smoke."""
from __future__ import annotations

import importlib.metadata
import importlib.util
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

from . import BASELINE, PROFILE, REGISTRY, ROOT, STUDY
from .contracts import read_yaml, study_issues
from .legacy import THREADS, evaluate_cell, load_profile, prepare_trusted
from .records import (event, identifier, new_directory, now, sha256, snapshot,
                      tree_hashes, verify_files, write_json)


def environment() -> dict:
    deps = {}
    for name in ("numpy", "scipy", "PyYAML", "numba", "pytest"):
        try:
            deps[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            deps[name] = None
    return {"python": sys.version, "executable": sys.executable, "dependencies": deps,
            "platform": platform.platform(), "processor": platform.processor(),
            "cpu_count": os.cpu_count(), "execution_threads": THREADS,
            "jit": "numba when available; original same-family warm-up",
            "cache_policy": "per-cell cache-free candidate; separate pristine cache",
            "os_isolation": "not_validated", "execution_policy": "reviewed builtins only"}


def doctor(study: Path = STUDY, profile: Path = PROFILE, registry: Path = REGISTRY,
           mode: str = "local") -> dict:
    problems = []
    try:
        pending = study_issues(read_yaml(study))
    except (OSError, ValueError, ImportError) as exc:
        pending = [str(exc)]
    try:
        resolved = load_profile(profile)
    except (OSError, ValueError, ImportError) as exc:
        problems.append(str(exc))
        resolved = None
    try:
        for method in read_yaml(registry)["methods"]:
            if method["implementation_kind"] == "archived_tree":
                problems.extend(verify_files({method["archive"]: method["archive_sha256"]}))
    except (OSError, ValueError, KeyError, TypeError, ImportError) as exc:
        problems.append(f"Registry: {exc}")
    env = environment()
    for name in ("numpy", "scipy", "PyYAML"):
        if env["dependencies"][name] is None:
            problems.append(f"Missing dependency: {name}")
    # Probe in a fresh process with the explicit repository source on sys.path.
    probe = subprocess.run([sys.executable, "-c",
        "import sys; sys.path.insert(0, sys.argv[1]); import pybte; print(pybte.__file__)",
        str(ROOT / "solver-python")], text=True, capture_output=True)
    actual = probe.stdout.strip()
    if probe.returncode or actual != str(ROOT / "solver-python/pybte/__init__.py"):
        problems.append("Repository pybte import failed or resolved to another installation")
    optional = {"mode": mode, "os_isolation": "not_validated"}
    if mode == "docker":
        optional.update(docker_executable=shutil.which("docker"), runtime_test="not_run",
                        container_claude_installation="not_demonstrated_by_legacy_Dockerfile")
    elif mode == "provider":
        optional.update(claude_executable=shutil.which("claude"), live_execution="blocked_by_C06")
    return {"local_smoke_ready": not problems, "environment": env, "pybte_import": actual,
            "profile_id": resolved["profile_id"] if resolved else None,
            "formal_execution_ready": False, "formal_blockers": pending,
            "formal_engine": "not_implemented_in_phase0", "input_errors": problems,
            "optional_mode": optional,
            "artifact_notice": "runs/ is ignored; use export for evidence worth retaining"}


def smoke(profile_path: Path, out: Path, method: str = "cis", provider: dict | None = None,
          timeout: float = 180) -> dict:
    if method not in ("cis", "krylov"):
        raise ValueError("Controlled smoke supports only reviewed builtin cis/krylov")
    if timeout <= 0:
        raise ValueError("timeout must be positive")
    profile = load_profile(profile_path)  # admission before any output or numerical work
    identifier(out.name)
    out = new_directory(out)
    started = now()
    event(out, "created", profile=profile["profile_id"])
    if provider:
        write_json(out / "provider.json", provider)
        event(out, "provider_finished", termination=provider["termination"])
    run = {"schema_version": 1, "run_id": out.name, "started_at": started,
           "purpose": "compatibility_smoke", "protocol_id": profile["protocol_id"],
           "counts_as_new_research_evidence": False,
           "scientific_claim_status": "not_evaluated_under_new_contracts",
           "environment": environment(), "provider": provider["provider"] if provider else "builtin",
           "new_protocol_evaluation": {"correctness_passed": None, "aggregate_score": None,
                                       "blocked_by": ["C01", "C02", "C03", "C04", "C05", "C06"]},
           "evaluations": [], "cost_usd": None, "export_status": "not_exported"}
    try:
        if provider and provider["termination"] != "completed":
            run.update(termination="provider_" + provider["termination"], legacy_gate=None)
        else:
            shutil.copyfile(profile_path, out / "resolved_config.yaml")
            candidate_dir = out / "candidates/submission_001"
            candidate = snapshot(candidate_dir, method, {"method": method})
            run["candidate"] = candidate
            event(out, "candidate_submitted", candidate_id=candidate["candidate_id"])
            trusted = prepare_trusted(out)
            frozen = tree_hashes(trusted)
            baseline = json.loads(BASELINE.read_text())
            write_json(out / "input_manifest.json", {
                "backend_id": "pybte_gray_callaway", "phase": "gray", "contract_status": "legacy_only",
                "source_revision": baseline["source_revision"], "verifier_revision": baseline["provenance"]["verifier_revision"],
                "profile_sha256": sha256(profile_path), "baseline_manifest_sha256": sha256(BASELINE),
                "trusted_files": frozen, "orchestration_files": tree_hashes(ROOT / "exploration"),
                "reference_certification": "per case: trusted/t01-square/reference/index.json; cis_iterations null means not cross-checked",
                "exposure": {"candidate_code": "candidates/submission_001/source", "method_hiding": False,
                             "agent_confirmation_or_answers": "not_exposed; no live provider"}})
            for i, name in enumerate(profile["cases"], 1):
                eval_id = f"eval_{i:03d}_{name}"
                folder = out / "evaluations" / eval_id
                folder.mkdir(parents=True, exist_ok=False)
                event(out, "execution_started", eval_id=eval_id, case=name)
                result, processes = evaluate_cell(trusted, candidate_dir / "source", name, method, folder, timeout)
                abnormal = next((p["termination"] for p in processes if p["termination"] != "completed"), None)
                status = abnormal or ("verifier_error" if "error" in result else "passed" if result["passed"] else "algorithm_failed")
                observation = {k: v for k, v in result.items() if k not in ("passed", "gates", "speedup", "scored")}
                write_json(folder / "observations.json", observation)
                evaluation = {"schema_version": 1, "eval_id": eval_id, "run_id": run["run_id"],
                              "candidate_id": candidate["candidate_id"], "case": name,
                              "protocol_id": profile["protocol_id"], "termination": status,
                              "legacy_gate": result.get("passed", False), "new_correctness_passed": None,
                              "new_aggregate_score": None, "raw_result": "raw_legacy.json",
                              "verifier_sha256": sha256(trusted / "verifier/evaluate.py"),
                              "gate_code_sha256": sha256(trusted / "verifier/checks.py"),
                              "reference_sha256": sha256(trusted / "reference" / f"{name}.npz"),
                              "artifacts": tree_hashes(folder)}
                # Work/cache trees are reproducible execution scratch, not evidence inputs.
                evaluation["artifacts"] = {k: v for k, v in evaluation["artifacts"].items()
                                           if not any("cache" in p or p == "candidate_work" for p in Path(k).parts)}
                write_json(folder / "evaluation.json", evaluation)
                run["evaluations"].append({"eval_id": eval_id, "path": str(folder.relative_to(out)),
                                           "case": name, "termination": status, "legacy_gate": evaluation["legacy_gate"]})
                event(out, "evaluation_finished", eval_id=eval_id, termination=status)
            if tree_hashes(candidate_dir / "source") != candidate["files"] or tree_hashes(trusted) != frozen:
                raise ValueError("Frozen candidate or trusted inputs changed during execution")
            original_issues = verify_files({p: h for p, h in baseline["files"].items()
                                           if p.startswith(tuple(baseline["protected_prefixes"]))})
            if original_issues:
                raise ValueError("Protected source changed: " + "; ".join(original_issues))
            run["legacy_gate"] = all(e["legacy_gate"] for e in run["evaluations"])
            failures = [e["termination"] for e in run["evaluations"] if e["termination"] != "passed"]
            run["termination"] = failures[0] if failures else "completed"
            run["integrity"] = "candidate, trusted copy and protected repository assets unchanged"
    except Exception as exc:
        run.update(termination="infrastructure_error", legacy_gate=False,
                   error=f"{type(exc).__name__}: {exc}")
        event(out, "failure", error=run["error"])
    run["finished_at"] = now()
    write_json(out / "run.json", run)
    lines = [f"# {run['run_id']}", "", f"Termination: {run['termination']}",
             f"Legacy gate: {run.get('legacy_gate')}", "",
             "Compatibility smoke only. New-protocol correctness and score are null; C01–C06 pending.",
             "OS sandbox: not_validated. Only reviewed builtins/mock executed.", "",
             "Temperature/heat flux are element integrals; vdf axes: direction, element, nodal DOF.",
             "The old gate checks temperature and transport residual; heat flux is recorded without a new gate.", "",
             "This run is in an ignored directory. Export evidence before relying on its preservation."]
    if run.get("error"):
        lines.extend(["", run["error"]])
    (out / "report.md").write_text("\n".join(lines) + "\n")
    event(out, "archived", termination=run["termination"])
    return run
