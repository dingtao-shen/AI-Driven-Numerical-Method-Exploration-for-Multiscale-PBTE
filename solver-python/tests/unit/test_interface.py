"""The public interface downstream tooling builds on.

These are contract tests, not behaviour tests.  They exist so that a rename or
a type change in `RunRecord` fails here rather than silently in a downstream
benchmark verifier months later.
"""
from __future__ import annotations

import inspect
import json

import numpy as np
import pytest

import pybte
from pybte import Case, RunRecord, Solver

from ..conftest import CASES


def test_top_level_names():
    """The three public names must be importable from the package root."""
    for name in ("Case", "Solver", "RunRecord"):
        assert hasattr(pybte, name)
    assert pybte.__version__


#: (attribute, expected type) of the frozen run-record interface
CONTRACT = [
    ("iterations", int),
    ("converged", bool),
    ("residual_history", np.ndarray),
    ("residual_true", (np.ndarray, type(None))),
    ("temp", np.ndarray),
    ("qx", np.ndarray),
    ("qy", np.ndarray),
    ("temp_dofs", np.ndarray),
    ("sweep_count", int),
    ("factorisation_count", int),
    ("wall_clock", float),
    ("config_hash", str),
]


@pytest.fixture(scope="module")
def record():
    c = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    c.flow.tau_r = 1.0
    c.velmesh.npole, c.velmesh.nazim, c.dg.deg = 6, 8, 1
    c.iteration.true_residual = True
    c.output.field = c.output.run_record = c.output.runtime_log = False
    return Solver(c).run()


@pytest.mark.parametrize("name,typ", CONTRACT)
def test_run_record_field(record, name, typ):
    assert hasattr(record, name), f"RunRecord lost {name}"
    assert isinstance(getattr(record, name), typ), name


def test_run_record_shapes(record):
    n_iter = record.iterations
    assert record.residual_history.shape == (n_iter,)
    assert record.residual_true.shape == (n_iter,)
    n_tris = record.n_tris
    assert record.temp.shape == (n_tris,)
    assert record.qx.shape == (n_tris,)
    assert record.qy.shape == (n_tris,)
    ndof = (record.deg + 1) * (record.deg + 2) // 2
    assert record.temp_dofs.shape == (ndof, n_tris)


def test_run_record_is_self_describing(record, tmp_path):
    """A verifier must never have to parse stdout, so everything needed to
    interpret the numbers has to be in the record itself."""
    p = record.to_json(tmp_path / "r.json")
    d = json.loads(p.read_text())
    for key in ("scheme", "tau_r", "tau_n", "tau_c", "deg", "npole", "nazim",
                "n_tris", "n_faces", "mesh_file", "tol", "tmax",
                "environment", "diagnostics", "restart",
                "final_residual", "mass"):
        assert key in d, f"run record does not carry {key}"
    assert d["environment"]["numpy"]
    assert d["scheme"] in ("CIS", "GSIS")
    assert len(d["config_hash"]) == 64


def test_sweep_count_is_hardware_independent(record):
    """sweep_count is the hardware-independent work metric; it must
    count sweeps, not seconds."""
    assert record.sweep_count == record.iterations


def test_solver_signature_is_stable():
    """`Solver(case)` must work with no other required argument."""
    sig = inspect.signature(Solver.__init__)
    required = [n for n, p in sig.parameters.items()
                if n != "self" and p.default is inspect.Parameter.empty]
    assert required == ["case"]
    assert inspect.signature(Solver.run).parameters["callback"].default is None


def test_case_from_yaml_is_a_classmethod():
    assert inspect.ismethod(Case.from_yaml)
