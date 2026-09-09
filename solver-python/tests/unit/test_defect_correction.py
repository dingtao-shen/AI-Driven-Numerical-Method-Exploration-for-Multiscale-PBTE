"""The experimental fixed-point repair -- docs/DEFECT_CORRECTION.md.

The property under test is the one the whole thing exists for: with the
defect on, GSIS lands on *CIS's* discrete fixed point rather than the
displaced one of the damped blend.  ``defect_omega = 0`` must leave the shipped
scheme untouched to the last bit.
"""
from __future__ import annotations

import numpy as np
import pytest

from pybte import Case, Solver

from .conftest_helpers import CASES, MESH_COARSE  # noqa: F401  (see below)


def _case(accflag, *, omega=0.0, every=1, anderson=0, tau_r=1.0,
          tol=1e-12, tmax=20000, deg=1, mesh="A1_Nx6_Ny6"):
    c = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    c.scheme.accflag = accflag
    c.scheme.defect_omega = omega
    c.scheme.defect_every = every
    c.scheme.defect_anderson = anderson
    c.flow.tau_r = tau_r
    c.dg.deg = deg
    c.mesh.file = f"../meshes/{mesh}.msh"
    c.iteration.tol = tol
    c.iteration.tmax = tmax
    c.output.field = c.output.run_record = c.output.runtime_log = False
    return c


def _run(**kw):
    s = Solver(_case(**kw))
    r = s.run()
    return s, r


def test_omega_zero_is_the_published_scheme():
    """Every defect operation is a no-op at omega = 0, to the last bit."""
    _, a = _run(accflag=1)
    _, b = _run(accflag=1, omega=0.0, every=0, anderson=8)
    assert a.iterations == b.iterations
    assert np.array_equal(a.temp, b.temp)
    assert a.mass == b.mass


@pytest.mark.slow
def test_defect_recovers_the_cis_fixed_point():
    """The point of the whole exercise.

    Published GSIS sits ~1e-2 away from CIS in the temperature field and
    leaves a transport residual of ~1e-3.  With the defect on, both collapse
    to round-off.
    """
    sc, rc = _run(accflag=0)
    sp, rp = _run(accflag=1)                                    # published
    sd, rd = _run(accflag=1, omega=1.0, every=0, anderson=8)    # repaired

    scale = np.abs(rc.temp).max()
    gap_pub = np.abs(rc.temp - rp.temp).max() / scale
    gap_fix = np.abs(rc.temp - rd.temp).max() / scale

    assert gap_pub > 1e-3, "the displacement should be there without the fix"
    assert gap_fix < 1e-8, f"defect correction left a gap of {gap_fix:.2e}"
    assert gap_fix < gap_pub / 1e4

    # and it is the *kinetic* fixed point, certified without reference to CIS
    assert sd.ctx.true_residual(sd.mom, sd.vdf) < 1e-10
    assert sp.ctx.true_residual(sp.mom, sp.vdf) > 1e-6


@pytest.mark.slow
def test_anderson_beats_plain_relaxation_on_the_outer_loop():
    """The outer map d -> o(d) is affine, so Anderson applies exactly."""
    sp, rp = _run(accflag=1, omega=1.0, every=0, anderson=0)
    sa, ra = _run(accflag=1, omega=1.0, every=0, anderson=8)
    assert len(sa.acc.defect_steps) < len(sp.acc.defect_steps)
    assert ra.iterations < rp.iterations


def test_config_validation():
    for bad in ({"defect_omega": -0.1}, {"defect_omega": 1.5},
                {"defect_every": -1}, {"defect_anderson": -1}):
        with pytest.raises(ValueError):
            from pybte.config import Scheme
            Scheme(**bad)


def test_diagnostics_report_the_defect():
    s, r = _run(accflag=1, omega=1.0, every=0, anderson=4)
    d = r.diagnostics
    assert d["defect_omega"] == 1.0
    assert d["defect_every"] == 0
    assert d["defect_anderson"] == 4
    assert d["defect_outer"] >= 1
    assert d["defect_norm"] > 0.0
