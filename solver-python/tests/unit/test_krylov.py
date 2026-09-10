"""``scheme.method: krylov`` -- same fixed point as source iteration, every wall type.

Cheap Knudsen numbers so the whole file runs in seconds; the property under
test does not depend on stiffness.
"""
from __future__ import annotations

import numpy as np
import pytest

from pybte import Case, Solver

from .conftest_helpers import CASES

CASE_FOR = {
    "isothermal": "cavity_tauR1e-3_cis.yaml",
    "adiabatic": "cavity_adiabatic_gsis.yaml",
    "periodic": "channel_periodic_gsis.yaml",
}


def _run(kind, method, tol, tmax=20000):
    c = Case.from_yaml(CASES / CASE_FOR[kind])
    c.scheme.accflag = 0
    c.scheme.method = method
    c.flow.tau_r = 1e-1
    c.flow.tau_n = 1e5
    c.dg.deg = 1
    c.velmesh.npole = 6
    c.velmesh.nazim = 8
    c.iteration.tol = tol
    c.iteration.tmax = tmax
    c.output.field = c.output.run_record = c.output.runtime_log = False
    s = Solver(c)
    return s, s.run()


@pytest.mark.parametrize("kind", ["isothermal", "adiabatic", "periodic"])
def test_krylov_lands_on_source_iterations_fixed_point(kind):
    sc, rc = _run(kind, "source", 1e-12)
    sk, rk = _run(kind, "krylov", 1e-8)
    assert rc.converged and rk.converged
    assert rk.scheme == "KRYLOV"
    gap = np.abs(rk.temp - rc.temp).max() / np.abs(rc.temp).max()
    assert gap < 1e-8, f"{kind}: Krylov sits {gap:.2e} from source iteration's fixed point"
    assert rk.diagnostics["transport_residual"] < 1e-10
    assert rk.iterations < rc.iterations


def test_periodic_state_is_part_of_the_krylov_vector():
    sk, rk = _run("periodic", "krylov", 1e-8)
    assert rk.diagnostics["krylov_periodic_state"] > 0
    # and the context is left the way source iteration expects it
    assert sk.ctx.use_pbuf is False


def test_source_iteration_never_uses_the_buffer():
    s, r = _run("periodic", "source", 1e-6, tmax=50)
    assert s.ctx.use_pbuf is False
