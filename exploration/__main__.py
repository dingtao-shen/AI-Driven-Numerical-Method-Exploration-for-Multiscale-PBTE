"""Run from the repository root: python -m exploration --help."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import BASELINE, PROFILE, REGISTRY, ROOT, STUDY
from .contracts import read_yaml, require_formal, study_issues
from .records import archive_members, export_run, safe_extract, sha256, write_json


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("doctor")
    p.add_argument("--study", type=Path, default=STUDY)
    p.add_argument("--profile", type=Path, default=PROFILE)
    p.add_argument("--registry", type=Path, default=REGISTRY)
    p.add_argument("--mode", choices=("local", "docker", "provider"), default="local")
    p = sub.add_parser("catalog")
    p.add_argument("--registry", type=Path, default=REGISTRY)
    p = sub.add_parser("run")
    p.add_argument("--study", type=Path, default=STUDY)
    p.add_argument("--dry-run", action="store_true")
    for name in ("smoke", "mock"):
        p = sub.add_parser(name)
        p.add_argument("--profile", type=Path, default=PROFILE)
        p.add_argument("--out", type=Path, required=True)
        p.add_argument("--method", choices=("cis", "krylov"), default="cis")
        p.add_argument("--timeout", type=float, default=180)
        if name == "mock":
            from .providers import SCENARIOS
            p.add_argument("--scenario", choices=SCENARIOS, default="success")
    p = sub.add_parser("recover")
    p.add_argument("--method", required=True)
    p.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("export")
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    try:
        if args.command == "doctor":
            from .pipeline import doctor
            result = doctor(args.study, args.profile, args.registry, args.mode)
            code = 0 if result["local_smoke_ready"] else 2
        elif args.command == "catalog":
            result, code = read_yaml(args.registry), 0
        elif args.command == "run":
            study = read_yaml(args.study)
            result = {"study_id": study.get("study_id"), "dry_run": args.dry_run,
                      "formal_blockers": study_issues(study), "execution": "not_started",
                      "prompt_template": "studies/gray_exploration/prompts/explore.md",
                      "exposure": {"knowledge": "requires C06", "editable_paths": "requires C02/C06",
                                   "budget": "requires C06", "evaluation": "requires C03/C04",
                                   "confirmation_inputs": "withheld", "reference_answers": "withheld"}}
            if not args.dry_run:
                require_formal(study)
                raise ValueError("Contracts structurally valid; formal research engine is not implemented in phase0")
            code = 0
        elif args.command in ("smoke", "mock"):
            from .pipeline import smoke
            provider = None
            if args.command == "mock":
                from .providers import mock_provider
                # Check output collision before invoking even a mock process.
                if args.out.exists():
                    raise FileExistsError(args.out)
                provider = mock_provider(args.scenario, args.method)
            result = smoke(args.profile, args.out, args.method, provider, args.timeout)
            code = 0 if result["termination"] == "completed" else 1
        elif args.command == "recover":
            entries = [m for m in read_yaml(REGISTRY)["methods"] if m["id"] == args.method]
            if len(entries) != 1 or entries[0]["implementation_kind"] != "archived_tree":
                raise ValueError("Recovery requires a registered archived_tree method")
            item = entries[0]
            from .records import new_directory
            out = new_directory(args.out)
            names = safe_extract(ROOT / item["archive"], out / "recovered", item["archive_sha256"])
            entry = out / "recovered/environment/pybte/driver.py"
            if not entry.is_file():
                raise ValueError("Archive does not contain the recorded driver entry point")
            result = {"method": args.method, "archive_sha256": item["archive_sha256"],
                      "members": len(names), "entrypoint": str(entry), "code_executed": False,
                      "numerical_replay": "not_run"}
            write_json(out / "recovery.json", result)
            code = 0
        else:
            manifest_path = args.out.with_suffix(args.out.suffix + ".manifest.json")
            if manifest_path.exists():
                raise FileExistsError(manifest_path)
            result = export_run(args.run.resolve(), args.out.resolve())
            write_json(manifest_path, result)
            code = 0
        print(json.dumps(result, indent=2, allow_nan=False))
        return code
    except (OSError, ValueError, KeyError, TypeError, ImportError) as exc:
        print(f"exploration: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
