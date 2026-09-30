import copy
import io
import json
import tarfile
from pathlib import Path

import pytest
import yaml

from exploration import BASELINE, PROFILE, ROOT, STUDY
from exploration.__main__ import main
from exploration.contracts import REQUIRED, digest_value, read_yaml, require_formal, study_issues
from exploration.legacy import load_profile
from exploration.records import (archive_members, export_run, identifier, new_directory,
                                 safe_extract, sha256, snapshot, tree_hashes, verify_files, write_json)


def approved_fixture():
    """Synthetic schema-only fixture; never a scientific approval or running profile."""
    study = read_yaml(STUDY)
    study.update(status="approved", formal_execution_enabled=True)
    for cid, fields in REQUIRED.items():
        value = {key: "test-only-structure" for key in fields}
        study["contracts"][cid] = {"status": "approved", "value": value,
            "approval": {"approved_by": "test-only", "approved_at": "2026-09-21T00:00:00Z",
                         "decision_ref": "test-fixture", "value_sha256": digest_value(value)}}
    return study


def test_pending_guard_and_no_boolean_bypass():
    study = read_yaml(STUDY)
    study["formal_execution_enabled"] = True
    with pytest.raises(ValueError, match="C01.*C06"):
        require_formal(study)
    assert study["non_gray_contract"]["status"] == "pending_owner"


@pytest.mark.parametrize("cid", ["C03", "C04", "C06"])
def test_missing_contracts_and_content_bound_approval(cid):
    study = approved_fixture()
    assert study_issues(study) == []  # C07 does not participate
    del study["contracts"][cid]
    assert any(cid in s for s in study_issues(study))
    study = approved_fixture()
    study["contracts"][cid]["value"][REQUIRED[cid][0]] = "different"
    assert f"{cid}: approval content hash mismatch" in study_issues(study)
    study["contracts"][cid]["approval"] = None
    assert f"{cid}: missing owner approval record" in study_issues(study)


@pytest.mark.parametrize("invalid", [None, {}, 0, "TODO", {"thresholds": None}])
def test_placeholder_contracts_are_not_approved(invalid):
    study = approved_fixture()
    study["contracts"]["C03"]["value"] = invalid
    assert any("C03" in s for s in study_issues(study))


def test_dry_run_never_executes_solver_or_provider(monkeypatch, capsys):
    import exploration.pipeline
    import exploration.providers
    def forbidden(*args, **kwargs):
        pytest.fail("dry-run executed work")
    monkeypatch.setattr(exploration.pipeline, "smoke", forbidden)
    monkeypatch.setattr(exploration.providers, "mock_provider", forbidden)
    before = sha256(STUDY)
    assert main(["run", "--dry-run"]) == 0
    assert '"reference_answers": "withheld"' in capsys.readouterr().out
    assert main(["run"]) == 2
    assert sha256(STUDY) == before


def test_bad_profile_missing_input_and_hash(tmp_path):
    profile = read_yaml(PROFILE)
    profile["profile_id"] = "unknown"
    p = tmp_path / "profile.yaml"
    p.write_text(yaml.safe_dump(profile))
    with pytest.raises(ValueError, match="profile_id"):
        load_profile(p)
    profile = read_yaml(PROFILE)
    profile["baseline_manifest_sha256"] = "0" * 64
    p.write_text(yaml.safe_dump(profile))
    with pytest.raises(ValueError, match="SHA256"):
        load_profile(p)
    assert verify_files({"missing": "0"}, tmp_path)
    assert verify_files({"profile.yaml": "0"}, tmp_path)


def test_exclusive_run_eval_ids_and_protected_outputs(tmp_path):
    out = new_directory(tmp_path / "run-001")
    with pytest.raises(FileExistsError):
        new_directory(out)
    write_json(out / "evaluation.json", {"test": 1})
    with pytest.raises(FileExistsError):
        write_json(out / "evaluation.json", {"test": 2})
    with pytest.raises(ValueError):
        identifier("../escape")
    with pytest.raises(ValueError, match="runs/"):
        new_directory(ROOT / "experiments/results/forbidden")


def test_snapshot_excludes_credentials_and_is_content_identified(tmp_path, monkeypatch):
    fake_home = tmp_path / "personal"
    fake_home.mkdir()
    (fake_home / ".env").write_text("TOKEN=not-a-real-secret")
    monkeypatch.setenv("HOME", str(fake_home))
    a = snapshot(tmp_path / "a", "cis", {"method": "cis"})
    b = snapshot(tmp_path / "b", "cis", {"method": "cis"})
    assert a["candidate_id"] == b["candidate_id"]
    assert all(p.startswith("pybte/") for p in a["files"])
    assert not any(".env" in p for p in a["files"])
    c = snapshot(tmp_path / "c", "krylov", {"method": "krylov"})
    assert c["candidate_id"] != a["candidate_id"]


@pytest.mark.parametrize("name,kind", [("../escape", "file"), ("/absolute", "file"),
                                      ("environment/link", "sym"), ("environment/hard", "hard")])
def test_unsafe_archive_rejected_before_extraction(tmp_path, name, kind):
    p = tmp_path / "candidate.tar.gz"
    with tarfile.open(p, "w:gz") as tf:
        item = tarfile.TarInfo(name)
        if kind in ("sym", "hard"):
            item.type = tarfile.SYMTYPE if kind == "sym" else tarfile.LNKTYPE
            item.linkname = "../../escape"
        tf.addfile(item)
    with pytest.raises(ValueError):
        safe_extract(p, tmp_path / "extracted", sha256(p))
    assert not (tmp_path / "extracted").exists()


def test_real_archive_recovery_and_wrong_hash(tmp_path):
    archive = ROOT / "experiments/results/t01-square-F3/trial_05_environment.tar.gz"
    expected = json.loads(BASELINE.read_text())["files"][str(archive.relative_to(ROOT))]
    with pytest.raises(ValueError, match="SHA256"):
        safe_extract(archive, tmp_path / "wrong", "0" * 64)
    names = safe_extract(archive, tmp_path / "restored", expected)
    assert "environment/pybte/driver.py" in names
    assert (tmp_path / "restored/environment/pybte/driver.py").is_file()


def test_export_contains_snapshot_and_can_be_recovered(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    write_json(run / "run.json", {"run_id": "example"})
    (run / ".env").write_text("DO_NOT_EXPORT=example")
    snapshot(run / "candidates/submission_001", "cis", {"method": "cis"})
    archive = tmp_path / "evidence.tar.gz"
    manifest = export_run(run, archive)
    assert ".env" not in manifest["files"]
    safe_extract(archive, tmp_path / "restore", manifest["sha256"])
    assert verify_files(manifest["files"], tmp_path / "restore") == []
    with pytest.raises(FileExistsError):
        export_run(run, archive)
    with pytest.raises(ValueError, match="protected"):
        export_run(run, ROOT / "experiments/results/forbidden.tar.gz")


def test_trusted_copy_uses_only_baseline_allowlist(tmp_path):
    from exploration.legacy import prepare_trusted
    baseline = json.loads(BASELINE.read_text())
    trusted = prepare_trusted(tmp_path)
    expected = {str(Path(p).relative_to("tasks/t01-square")): h
                for p, h in baseline["files"].items() if p.startswith("tasks/t01-square/")}
    assert tree_hashes(trusted) == expected
