"""§6.3 integration level: full runs against the Fortran golden logs.

Iteration counts must match *exactly*; the converged field to the per-gate
tolerance.  Where the Fortran log exists these compare against it; the rest
are self-consistency properties that hold with no Fortran at all.
"""
from __future__ import annotations

import numpy as np
import pytest

from pybte import Case, Solver

from ..conftest import CASES, GOLDEN, has_golden

pytestmark = pytest.mark.slow


def _case(tau_r=1e-3, accflag=0, deg=3, tol=1e-8, tmax=8_000_000, mesh=None):
    c = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    c.scheme.accflag = accflag
    c.flow.tau_r = tau_r
    c.dg.deg = deg
    c.iteration.tol = tol
    c.iteration.tmax = tmax
    if mesh:
        c.mesh.file = f"../meshes/{mesh}.msh"
    c.output.field = c.output.run_record = c.output.runtime_log = False
    return c


def _golden(name):
    if not has_golden(name):
        pytest.skip(f"no Fortran golden run {name}")
    return np.atleast_2d(np.loadtxt(GOLDEN / name / "residual_history.txt"))


@pytest.mark.fortran
@pytest.mark.parametrize("tau_r,name", [
    (1.0, "cis_tauR1.0_deg3"),
    (1e-1, "cis_tauR1e-1_deg3"),
])
def test_cis_residual_history_is_bit_identical(tau_r, name):
    """Gate 3 asks for rtol=1e-10 over 50 iterations; we get bit equality over
    the whole run, which is a far stronger statement about the port."""
    ref = _golden(name)
    rec = Solver(_case(tau_r=tau_r)).run()
    assert rec.iterations == int(ref[-1, 0])
    assert np.array_equal(rec.residual_history, ref[:, 1])
    assert np.array_equal(rec.mass_history, ref[:, 2])


@pytest.mark.fortran
def test_gsis_matches_fortran_iteration_count():
    ref = _golden("gsis_tauR1e-3_deg3")
    rec = Solver(_case(tau_r=1e-3, accflag=1)).run()
    assert rec.iterations == int(ref[-1, 0])
    n = min(len(ref), rec.iterations)
    rel = np.abs(rec.residual_history[:n] - ref[:n, 1]) / np.abs(ref[:n, 1])
    assert rel.max() < 1e-6, f"residual history drifted by {rel.max():.2e}"


# ---------------------------------------------------------------------------
# properties that need no Fortran
# ---------------------------------------------------------------------------
def test_cis_reaches_the_discrete_fixed_point():
    """The stopping criterion is on the *iterate* difference, but the solution
    it stops at genuinely solves the discrete transport system."""
    c = _case(tau_r=1e-1, tol=1e-12)
    s = Solver(c)
    s.run()
    assert s.ctx.true_residual(s.mom, s.vdf) < 1e-10


def test_gsis_is_far_cheaper_than_cis():
    """The phenomenon the whole project is about: CIS iteration count grows by
    orders of magnitude as Kn falls, GSIS stays roughly flat."""
    counts = {}
    for tau_r in (1e-1, 1e-2):
        counts[("cis", tau_r)] = Solver(_case(tau_r=tau_r)).run().iterations
        counts[("gsis", tau_r)] = Solver(_case(tau_r=tau_r, accflag=1)).run().iterations
    assert counts[("cis", 1e-2)] > 10 * counts[("cis", 1e-1)]
    assert counts[("gsis", 1e-2)] < 3 * counts[("gsis", 1e-1)]
    assert counts[("gsis", 1e-2)] < counts[("cis", 1e-2)] / 100


def test_cis_mass_approaches_the_diffusion_limit():
    """``int T dA`` over the unit square is exactly 1/4 for the Laplace
    solution, and the kinetic solution approaches it from below as Kn falls.

    CIS only.  GSIS is *not* monotone here -- its mass error goes 3.5e-6 at
    ``tau_R = 1e-1`` but 3.8e-5 at ``1e-2`` -- because it does not converge to
    the kinetic fixed point at all (docs/FORTRAN_ISSUES.md #5).  Asserting
    monotonicity for GSIS would be asserting something false.
    """
    prev = None
    for tau_r in (1.0, 1e-1, 1e-2):
        rec = Solver(_case(tau_r=tau_r, accflag=0)).run()
        assert rec.converged
        assert 0.24 < rec.mass < 0.25
        if prev is not None:
            assert abs(rec.mass - 0.25) < abs(prev - 0.25)
        prev = rec.mass


def test_gsis_mass_is_close_but_not_convergent():
    """GSIS stays within 2e-4 of the diffusion limit at every Kn, which is
    useful as a sanity band -- but see the previous test for why it is stated
    as a bound rather than as convergence."""
    for tau_r in (1.0, 1e-1, 1e-2):
        rec = Solver(_case(tau_r=tau_r, accflag=1)).run()
        assert rec.converged
        assert abs(rec.mass - 0.25) < 2e-4


def test_run_record_is_complete():
    """§8: a verifier must never need to parse stdout."""
    rec = Solver(_case(tau_r=1.0, accflag=1)).run()
    for f in ("iterations", "converged", "residual_history", "temp", "qx", "qy",
              "temp_dofs", "sweep_count", "factorisation_count", "wall_clock",
              "config_hash"):
        assert getattr(rec, f) is not None, f
    assert rec.sweep_count == rec.iterations
    assert rec.scheme == "GSIS" and rec.acc_variant == "A"
    assert len(rec.config_hash) == 64
    d = rec.to_dict()
    import json
    json.dumps(d)                       # must be JSON-serialisable


def test_deg_convergence_towards_analytic():
    """Raising the polynomial order moves the field towards the Fourier limit
    at small Kn -- an independent check that does not involve the Fortran."""
    from pybte.analytic import fourier_temperature

    errs = []
    for deg in (1, 2, 3):
        c = _case(tau_r=1e-2, accflag=1, deg=deg)
        s = Solver(c)
        rec = s.run()
        cen = s.mesh.nodes[s.mesh.tri_nodes].mean(axis=1)
        tavg = rec.temp / s.mesh.tri_area
        ref = fourier_temperature(cen[:, 0], cen[:, 1])
        interior = np.minimum.reduce(
            [cen[:, 0], 1 - cen[:, 0], cen[:, 1], 1 - cen[:, 1]]) > 0.1
        errs.append(float(np.sqrt(np.mean((tavg - ref)[interior] ** 2))))
    assert errs[2] < errs[0], f"no improvement with DEG: {errs}"
    assert errs[2] < 0.02
