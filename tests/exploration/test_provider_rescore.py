import importlib.util
import json
import shutil
import sys

import pytest

from exploration import ROOT
from exploration.providers import mock_provider
from exploration.records import tree_hashes


@pytest.mark.parametrize("scenario,expected", [("success", "completed"), ("timeout", "timeout"),
    ("error", "execution_error"), ("empty_transcript", "empty_transcript"),
    ("missing_final", "missing_final_event"), ("invalid_output", "invalid_output")])
def test_mock_lifecycle(scenario, expected):
    result = mock_provider(scenario)
    assert result["termination"] == expected
    assert result["model_id"] is None and result["cost_usd"] is None
    if expected == "completed":
        assert result["events"][0] == {"type": "candidate", "method": "cis"}


def rescore_module():
    spec = importlib.util.spec_from_file_location("rescore", ROOT / "tools/rescore_rollouts.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_rescore_append_and_original_hashes(tmp_path, monkeypatch):
    module = rescore_module()
    source = ROOT / "experiments/results/t01-square-F3"
    before = tree_hashes(source)
    calls = []
    def fake_grade(task, env, output, timeout):
        assert (env / "pybte/driver.py").is_file()
        calls.append(output)
        return {"gate": True, "score": 1.0, "families": {}, "termination": "passed"}
    monkeypatch.setattr(module, "grade", fake_grade)
    args = [str(source), "--task", str(ROOT / "tasks/t01-square"),
            "--sandbox-root", str(tmp_path / "scratch"), "--trials", "5"]
    assert module.main(args + ["--out", str(tmp_path / "eval-001")]) == 0
    assert module.main(args + ["--out", str(tmp_path / "eval-002")]) == 0
    assert module.main(args + ["--out", str(tmp_path / "eval-001")]) == 2
    assert len(calls) == 2
    assert tree_hashes(source) == before
    meta = json.loads((tmp_path / "eval-002/evaluation_manifest.json").read_text())
    assert meta["source_unchanged"] is True
    assert meta["counts_as_new_research_evidence"] is False
    assert module.main(args + ["--in-place-copy"]) == 2


def test_rescore_missing_archive_duplicates_and_copy_links(tmp_path):
    module = rescore_module()
    args = [str(ROOT / "experiments/results/t01-square-F3"), "--task", str(ROOT / "tasks/t01-square"),
            "--sandbox-root", str(tmp_path / "scratch"), "--out", str(tmp_path / "out")]
    assert module.main(args + ["--trials", "99"]) == 2
    assert module.main(args + ["--trials", "1", "1"]) == 2
    assert not (tmp_path / "out").exists()
    copy = tmp_path / "copy"
    copy.mkdir()
    (copy / "linked").symlink_to(ROOT / "experiments/results/t01-square-F3/summary.json")
    with pytest.raises(ValueError, match="links"):
        module.check_writable_copy(copy)


def test_legacy_rollout_requires_fresh_explicit_output_without_model_calls(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("rollouts", ROOT / "tools/run_rollouts.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    def forbidden(*args, **kwargs):
        pytest.fail("guard called a live provider")
    monkeypatch.setattr(module, "run_agent_dir", forbidden)
    args = [str(ROOT / "tasks/t01-square")]
    assert module.main(args) == 2
    assert module.main(args + ["--out", str(ROOT / "experiments/results/t01-square-F3")]) == 2
    assert module.main(args + ["--out", str(tmp_path)]) == 2
