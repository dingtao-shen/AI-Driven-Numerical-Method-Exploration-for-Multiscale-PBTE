import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

from exploration import PROFILE, ROOT, STUDY
from exploration.contracts import read_yaml
from exploration.legacy import THREADS, run_stage, validate_solution
from exploration.pipeline import smoke
from exploration.providers import mock_provider
from exploration.records import sha256


@pytest.fixture(scope="module")
def numerical_run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("numerical")
    before = sha256(STUDY)
    result = smoke(PROFILE, tmp / "actual-smoke", "krylov", mock_provider("success", "krylov"))
    assert sha256(STUDY) == before
    assert result["termination"] == "completed", result
    return tmp / "actual-smoke", result


def test_real_mock_submission_solve_evaluate(numerical_run):
    out, result = numerical_run
    assert result["legacy_gate"] is True
    assert result["new_protocol_evaluation"]["correctness_passed"] is None
    assert len(result["evaluations"]) == 3
    assert [e["case"] for e in result["evaluations"]] == ["F1_1_1", "F2_1_1", "F3_1_1"]
    for e in result["evaluations"]:
        path = out / e["path"]
        raw = json.loads((path / "raw_legacy.json").read_text())
        assert raw["warm_up"] == "ok"
        assert raw["transport_residual"] < 1e-7
        assert raw["field_gap_worst"] < 1e-5
        details = json.loads((path / f"{e['case']}.details.json").read_text())
        verifier = json.loads((path / f"{e['case']}.verifier_import.json").read_text())
        assert "candidate_work/pybte" in details["pybte_import"]
        assert "verifier/pristine/pybte" in verifier["pybte_import"]
        with np.load(path / f"{e['case']}.npz") as data:
            assert data["vdf"].shape == (200, 200, 6)
            assert data["temp"].shape == data["qx"].shape == data["qy"].shape == (200,)


@pytest.mark.parametrize("mode", ["shape", "temp_shape", "nan", "inf", "forged_convergence"])
def test_independent_verifier_rejects_corrupt_or_forged_solution(numerical_run, tmp_path, mode):
    out, result = numerical_run
    cell = out / result["evaluations"][0]["path"]
    case = cell / "F1_1_1.yaml"
    with np.load(cell / "F1_1_1.npz") as original:
        arrays = {k: original[k] for k in ("vdf", "temp")}
    if mode == "shape":
        arrays["vdf"] = arrays["vdf"][:, :-1, :]
    elif mode == "temp_shape":
        arrays["temp"] = arrays["temp"][:-1]
    elif mode in ("nan", "inf"):
        arrays["vdf"][0, 0, 0] = float(mode)
    else:
        arrays["vdf"].fill(0)
        arrays["converged"] = True
        arrays["true_residual"] = 0.0  # must not be trusted
    dump = tmp_path / "bad.npz"
    np.savez(dump, **arrays)
    if mode in ("nan", "inf"):
        with pytest.raises(ValueError, match="non-finite"):
            validate_solution(dump)
        return
    trusted = out / "trusted/t01-square"
    command = [sys.executable, str(ROOT / "exploration/worker.py"), "evaluate",
               str(trusted / "verifier/evaluate.py"), str(trusted / "verifier/pristine"),
               str(case), str(dump), str(trusted / "reference/F1_1_1.npz")]
    data, err, meta = run_stage(command, tmp_path, "independent", 120,
                                dict(os.environ, **THREADS, NUMBA_CACHE_DIR=str(tmp_path / "cache")))
    assert err is None
    if mode in ("shape", "temp_shape"):
        assert "shape" in data["error"]
    else:
        assert data["transport_residual"] > 1e-7
        assert data["field_gap_from_vdf"] > 1e-5


@pytest.mark.parametrize("scenario", ["timeout", "error", "empty_transcript", "missing_final", "invalid_output"])
def test_provider_failures_are_archived_without_solving(tmp_path, scenario):
    result = smoke(PROFILE, tmp_path / scenario, provider=mock_provider(scenario))
    assert result["termination"].startswith("provider_")
    assert result["legacy_gate"] is None
    assert result["evaluations"] == []
    assert (tmp_path / scenario / "run.json").exists()


def test_warmup_failure_never_scored(monkeypatch, tmp_path):
    from exploration import legacy
    def failed_stage(*args, **kwargs):
        return {"warm_up": "failed intentionally"}, None, {"termination": "completed"}
    monkeypatch.setattr(legacy, "run_stage", failed_stage)
    from exploration.records import snapshot
    source = tmp_path / "source"
    snapshot(source, "cis", {"method": "cis"})
    trusted = legacy.prepare_trusted(tmp_path)
    output = tmp_path / "evaluation"
    output.mkdir()
    result, _ = legacy.evaluate_cell(trusted, source / "source", "F1_1_1", "cis", output)
    assert result["passed"] is False
    assert "speedup" not in result and "sweep_equivalents" not in result


@pytest.mark.parametrize("program,expected", [("import time; time.sleep(5)", "timeout"),
    ("import sys; sys.exit(9)", "execution_error"), ("print('broken')", "invalid_output")])
def test_subprocess_failure_status_and_logs(tmp_path, program, expected):
    data, error, meta = run_stage([sys.executable, "-c", program], tmp_path, "failure", 0.1, os.environ.copy())
    assert data is None and error
    assert meta["termination"] == expected
    assert json.loads((tmp_path / "failure.process.json").read_text())["termination"] == expected
    assert (tmp_path / "failure.stdout.log").exists()
