"""Properties of a full run.

Numerical agreement with the Fortran reference is frozen in
``test_regression.py``; these are the physical and structural properties that
have to hold on their own.
"""
from __future__ import annotations

import numpy as np
import pytest

from pybte import Case, Solver

from ..conftest import CASES

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
    the kinetic fixed point at all (docs/LIMITATIONS.md #1).  Asserting
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


def test_gsis_mass_tracks_cis_not_the_diffusion_limit():
    """GSIS's `int T dA` must be judged against CIS at the same Kn, not
    against the diffusion limit.

    Comparing to 0.25 conflates two different things.  At `tau_R = 1` the
    kinetic answer genuinely is not 0.25 -- both schemes give 0.249687, and
    the 3.1e-4 shortfall is physics, not error.  What is worth bounding is how
    far GSIS lands from CIS, which is the fixed-point gap of
    docs/LIMITATIONS.md #1, and that shrinks as Kn falls.
    """
    gaps = {}
    for tau_r in (1.0, 1e-1, 1e-2):
        cis = Solver(_case(tau_r=tau_r, accflag=0)).run()
        gsis = Solver(_case(tau_r=tau_r, accflag=1)).run()
        assert cis.converged and gsis.converged
        gaps[tau_r] = abs(gsis.mass - cis.mass)
        assert gaps[tau_r] < 1e-3, f"tau_R={tau_r:g}: gap {gaps[tau_r]:.2e}"
    assert gaps[1e-2] < gaps[1e-1], f"gap should shrink with Kn: {gaps}"


def test_run_record_is_complete():
    """A verifier must never need to parse stdout."""
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
