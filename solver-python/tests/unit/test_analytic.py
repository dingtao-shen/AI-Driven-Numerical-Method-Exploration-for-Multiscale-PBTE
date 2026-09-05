"""The Fourier (diffusion-limit) reference solution."""
from __future__ import annotations

import numpy as np
import pytest

from pybte.analytic import (analytic_fields, fourier_flux, fourier_temperature,
                            l2_error)
from pybte.constants import PI


def test_boundary_values():
    """Three of the four walls are satisfied exactly term by term."""
    x = np.linspace(0.001, 0.999, 41)
    assert np.allclose(fourier_temperature(x, np.zeros_like(x)), 0.0, atol=1e-14)
    y = np.linspace(0.0, 0.999, 41)
    assert np.allclose(fourier_temperature(np.zeros_like(y), y), 0.0, atol=1e-14)
    assert np.allclose(fourier_temperature(np.ones_like(y), y), 0.0, atol=1e-13)


def test_hot_wall_is_one_up_to_gibbs():
    """On the hot wall the truncated series rings.

    The boundary datum is discontinuous at the two top corners, so a 200-term
    truncation overshoots by ~3% next to them and only settles to 1 in the
    middle.  This is a property of the reference solution the Fortran ships,
    not of the port -- it bounds how sharply any solver can be compared
    against it near y = 1.
    """
    x = np.linspace(0.2, 0.8, 31)
    t = fourier_temperature(x, np.ones_like(x))
    assert np.allclose(t, 1.0, atol=1.2e-2)
    assert np.abs(t - 1.0).mean() < 5e-3
    # ... and the overshoot really is there near the corner
    assert fourier_temperature(0.02595, 1.0) > 1.0


def test_mean_is_one_quarter():
    """By symmetry, superposing the four rotations of this problem gives the
    constant 1, so the mean of a single one is exactly 1/4."""
    n = 400
    t = (np.arange(n) + 0.5) / n
    X, Y = np.meshgrid(t, t, indexing="ij")
    assert float(fourier_temperature(X, Y).mean()) == pytest.approx(0.25, abs=1e-4)


def test_satisfies_laplace():
    """The series is harmonic in the interior.

    Checked by second-order differences at two step sizes: what is left is the
    ``h^2/12 * grad^4 T`` truncation of the stencil, so the residual has to
    fall by ~4 when ``h`` halves.  Asserting the *rate* rather than a fixed
    tolerance is what makes this a test of the series and not of the stencil.
    """
    x = np.array([0.31, 0.5, 0.77])
    y = np.array([0.29, 0.5, 0.61])
    X, Y = np.meshgrid(x, y, indexing="ij")

    def lap(h):
        return (fourier_temperature(X + h, Y) + fourier_temperature(X - h, Y)
                + fourier_temperature(X, Y + h) + fourier_temperature(X, Y - h)
                - 4 * fourier_temperature(X, Y)) / h ** 2

    coarse = np.abs(lap(4e-3)).max()
    fine = np.abs(lap(1e-3)).max()
    assert coarse < 1e-3
    assert fine < 1e-4
    assert coarse / fine == pytest.approx(16.0, rel=0.3)


def test_flux_is_fourier_law():
    """q = -(Cv Vg^2 tau/3) grad T, checked against a central difference."""
    cv, vg, tau = 1.3, 0.7, 2.5e-3
    kappa = cv * vg * vg * tau / 3.0
    h = 1e-5
    x, y = 0.37, 0.62
    dTdx = (fourier_temperature(x + h, y) - fourier_temperature(x - h, y)) / (2 * h)
    dTdy = (fourier_temperature(x, y + h) - fourier_temperature(x, y - h)) / (2 * h)
    qx, qy = fourier_flux(x, y, cv, vg, tau)
    assert float(qx) == pytest.approx(-kappa * dTdx, rel=1e-6)
    assert float(qy) == pytest.approx(-kappa * dTdy, rel=1e-6)


def test_fortran_compat_flux_is_half():
    """Out_Put_Result.f90 differentiates the series but forgets the surviving
    factor 2, so its Conduction_A.dat fluxes are half the correct value
    (docs/FORTRAN_ISSUES.md #3).  Keep the switch, default to correct."""
    a = fourier_flux(0.4, 0.6, 1.0, 1.0, 1e-3, fortran_compat=False)
    b = fourier_flux(0.4, 0.6, 1.0, 1.0, 1e-3, fortran_compat=True)
    assert float(a[0]) == pytest.approx(2.0 * float(b[0]), rel=1e-14)
    assert float(a[1]) == pytest.approx(2.0 * float(b[1]), rel=1e-14)


def test_l2_error_utility():
    t = (np.arange(30) + 0.5) / 30
    X, Y = np.meshgrid(t, t, indexing="ij")
    exact = fourier_temperature(X, Y)
    e = l2_error(X, Y, exact)
    assert e["l2_abs"] == pytest.approx(0.0, abs=1e-14)
    e2 = l2_error(X, Y, exact + 0.01)
    assert e2["l2_abs"] == pytest.approx(0.01, rel=1e-10)


def test_series_truncation_guard():
    with pytest.raises(ValueError):
        fourier_temperature(0.5, 0.5, n_terms=500)
