"""Exclusive creation, full hashes, allowlisted snapshots and safe archives."""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tarfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

from . import BASELINE, ROOT


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value) -> None:
    with path.open("x") as fh:
        json.dump(value, fh, indent=2, allow_nan=False)
        fh.write("\n")


def new_directory(path: Path) -> Path:
    path = path.resolve()
    # New outputs must never be nested inside source, fixtures or evidence.
    if path == ROOT or (ROOT in path.parents and path.relative_to(ROOT).parts[0] not in ("runs", "workspaces")):
        raise ValueError("Repository outputs must be under runs/ or workspaces/")
    path.mkdir(parents=True, exist_ok=False)
    return path


def identifier(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,100}", value):
        raise ValueError("Invalid run/eval/candidate ID")
    return value


def event(path: Path, kind: str, **data) -> None:
    with (path / "events.jsonl").open("a") as fh:
        fh.write(json.dumps({"time": now(), "event": kind, **data}, allow_nan=False) + "\n")


def tree_hashes(path: Path) -> dict:
    return {str(p.relative_to(path)): sha256(p) for p in sorted(path.rglob("*"))
            if p.is_file() and "__pycache__" not in p.parts}


def verify_files(files: dict, root: Path = ROOT) -> list[str]:
    issues = []
    for name, expected in files.items():
        p = root / name
        if not p.is_file() or p.is_symlink():
            issues.append(f"Missing input or symlink: {name}")
        elif sha256(p) != expected:
            issues.append(f"SHA256 mismatch: {name}")
    return issues


def snapshot(out: Path, method: str, config: dict) -> dict:
    source = out / "source"
    source.mkdir(parents=True, exist_ok=False)
    baseline = json.loads(BASELINE.read_text())
    # Explicit known source files only: never copy HOME, credentials, caches or logs.
    files = {p: h for p, h in baseline["files"].items() if p.startswith("solver-python/pybte/")}
    issues = verify_files(files)
    if issues:
        raise ValueError("; ".join(issues))
    for name in files:
        target = source / Path(name).relative_to("solver-python")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
        target.chmod(0o444)
    write_json(out / "method_config.json", config)
    code_hashes = tree_hashes(source)
    candidate_id = "candidate_" + hashlib.sha256(json.dumps(
        {"files": code_hashes, "config": config}, sort_keys=True).encode()).hexdigest()
    result = {"schema_version": 1, "candidate_id": candidate_id, "parent_method": method,
              "parent_candidate": None, "source_revision": baseline["source_revision"],
              "files": code_hashes, "config_sha256": sha256(out / "method_config.json"),
              "snapshot_policy": "allowlisted source; read-only files; hashes checked before and after execution"}
    write_json(out / "manifest.json", result)
    return result


def archive_members(archive: Path) -> list[str]:
    """Reject links, devices, traversal, duplicate files and oversized archives."""
    names, seen, total = [], set(), 0
    with tarfile.open(archive) as tf:
        for member in tf:
            p = PurePosixPath(member.name)
            if p.is_absolute() or ".." in p.parts or "\\" in member.name:
                raise ValueError(f"Unsafe archive path: {member.name}")
            if not (member.isfile() or member.isdir()):
                raise ValueError(f"Archive links/special files are not allowed: {member.name}")
            normalized = str(p)
            if normalized in seen:
                raise ValueError(f"Duplicate archive path: {member.name}")
            seen.add(normalized)
            total += member.size
            if len(seen) > 100000 or total > 2 * 1024**3:
                raise ValueError("Archive exceeds engineering recovery limits")
            names.append(normalized)
    return names


def safe_extract(archive: Path, destination: Path, expected_hash: str) -> list[str]:
    if sha256(archive) != expected_hash:
        raise ValueError("Archive SHA256 mismatch")
    names = archive_members(archive)
    destination.mkdir(parents=True, exist_ok=False)
    with tarfile.open(archive) as tf:
        tf.extractall(destination, filter="data")
    return names


def clean_copy(source: Path, target: Path) -> None:
    shutil.copytree(source, target, ignore=shutil.ignore_patterns(
        "__pycache__", "*.pyc", "*.nbc", "*.nbi", ".pytest_cache", "out"))
    for p in target.rglob("*"):
        if p.is_file():
            p.chmod(0o644)


def export_run(source: Path, target: Path) -> dict:
    """Export known run artifacts, excluding execution work trees and caches."""
    target = target.resolve()
    for protected in (ROOT / "tasks", ROOT / "experiments/results", ROOT / "solver-python"):
        if target == protected or protected in target.parents:
            raise ValueError("Export cannot write into protected source/evidence directories")
    run = json.loads((source / "run.json").read_text())
    prefixes = ("candidates", "evaluations", "trusted", "inputs")
    root_files = {"run.json", "resolved_config.yaml", "input_manifest.json", "events.jsonl",
                  "report.md", "provider.json"}
    files = []
    for p in sorted(source.rglob("*")):
        rel = p.relative_to(source)
        if p.is_file() and ((len(rel.parts) == 1 and p.name in root_files) or rel.parts[0] in prefixes):
            if any("cache" in part or part == "candidate_work" for part in rel.parts) or p.suffix in (".pyc", ".nbi", ".nbc"):
                continue
            if p.is_symlink() or p.name.startswith(".env") or p.suffix in (".pem", ".key"):
                raise ValueError(f"Unsafe export artifact: {rel}")
            files.append(p)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("xb") as fh, tarfile.open(fileobj=fh, mode="w:gz") as tf:
        for p in files:
            tf.add(p, arcname=str(p.relative_to(source)), recursive=False)
    return {"run_id": run["run_id"], "archive": str(target), "sha256": sha256(target),
            "files": {str(p.relative_to(source)): sha256(p) for p in files},
            "recovery": "Extract with safe_extract into a fresh directory; verify files against this manifest. Use the baseline revision and pinned legacy profile for replay; recorded timings are machine dependent."}
