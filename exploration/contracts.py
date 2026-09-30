"""Fail-closed structure checks. Scientific approval belongs to the owner."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

REQUIRED = {
    "C01": ("model", "kn_domain", "geometry", "boundaries", "forcing", "normalization"),
    "C02": ("discretization", "allowed_changes", "file_permissions"),
    "C03": ("residual", "field_metrics", "thresholds", "zero_scale", "reference_certification"),
    "C04": ("primary_metric", "cost_accounting", "cis_calibration", "censoring", "repetitions", "aggregation", "acceleration_criterion"),
    "C05": ("exploration_split", "confirmation_split", "resolution_checks", "seeds", "candidate_freeze", "controls"),
    "C06": ("knowledge", "editable_paths", "provider", "model", "scaffold", "budget", "network", "dependencies", "isolation", "abnormal_trials"),
}


def read_yaml(path: Path) -> dict:
    import yaml
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected a YAML mapping")
    return data


def digest_value(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def complete(value) -> bool:
    if value is None or value == "" or value == {} or value == []:
        return False
    if isinstance(value, str):
        return value.strip().lower() not in ("", "todo", "tbd", "pending", "pending_owner")
    if isinstance(value, dict):
        return all(complete(v) for v in value.values())
    if isinstance(value, list):
        return all(complete(v) for v in value)
    return isinstance(value, (bool, int, float))


def study_issues(study: dict) -> list[str]:
    issues = []
    if study.get("schema_version") != 1 or study.get("phase") != "gray":
        issues.append("schema_version must be 1 and phase must be gray")
    if not isinstance(study.get("study_id"), str) or not study["study_id"].strip():
        issues.append("study_id is required")
    if study.get("status") != "approved":
        issues.append("study status is not approved")
    if study.get("formal_execution_enabled") is not True:
        issues.append("formal_execution_enabled is not true")
    contracts = study.get("contracts")
    if not isinstance(contracts, dict):
        contracts = {}
    for cid, fields in REQUIRED.items():
        item = contracts.get(cid)
        if not isinstance(item, dict) or item.get("status") != "approved":
            issues.append(f"{cid}: pending_owner or missing approval")
            continue
        value = item.get("value")
        if not isinstance(value, dict) or not all(k in value and complete(value[k]) for k in fields):
            issues.append(f"{cid}: incomplete value; required fields: {', '.join(fields)}")
            continue
        approval = item.get("approval")
        if not isinstance(approval, dict) or not all(
            isinstance(approval.get(k), str) and complete(approval[k])
            for k in ("approved_by", "approved_at", "decision_ref", "value_sha256")
        ):
            issues.append(f"{cid}: missing owner approval record")
            continue
        try:
            if datetime.fromisoformat(approval["approved_at"].replace("Z", "+00:00")).tzinfo is None:
                raise ValueError("timezone required")
            if approval["value_sha256"] != digest_value(value):
                issues.append(f"{cid}: approval content hash mismatch")
        except (ValueError, TypeError):
            issues.append(f"{cid}: invalid approval timestamp or value")
    return issues


def require_formal(study: dict) -> None:
    issues = study_issues(study)
    if issues:
        raise ValueError("Formal execution blocked before solve/provider: " + "; ".join(issues))
