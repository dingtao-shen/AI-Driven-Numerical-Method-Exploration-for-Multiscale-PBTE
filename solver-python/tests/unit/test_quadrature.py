"""Quadrature rules and the angular mesh."""
from __future__ import annotations

import numpy as np
import pytest

from pybte.constants import PI
from pybte.quadrature import (gauss_legendre, leggauss_reference,
                              tri_quadrature)
from pybte.velocity import build_velocity_mesh


@pytest.mark.parametrize("n", [2, 3, 5, 10, 15, 20, 32])
def test_gauss_legendre_matches_numpy(n):
    """The Fortran's hand-rolled iteration lands on the right rule, but only
    to ~1e-13 absolute -- not to machine precision.

    The reason is in ``USD_Math.f90``: it forms

        Lp = N2*(L(:,N1) - y*L(:,N2)) / (1 - y*y)

    with ``N2 = n+1``, whereas the Legendre identity has ``n``.  The weight
    formula carries the same ``N2**2/N1**2`` factor, so the *weights* are
    consistent with whatever roots come out -- but the Newton step is scaled
    by ``n/(n+1)``, which turns quadratic convergence into linear convergence
    with ratio ``1/(n+1)``.  Combined with the step-based exit test
    (``|dy| < 1e-13``) the roots stop about ``1e-13/n`` short.

    Downstream this shows up as ``sum(DOMEGA)`` missing ``4 pi`` by ~2e-14
    relative.  It is reference behaviour, reproduced bit-for-bit, and is
    recorded in docs/LIMITATIONS.md #4 rather than silently fixed.
    """
    a, w = gauss_legendre(n, 0.0, PI)
    ra, rw = leggauss_reference(n, 0.0, PI)
    assert np.allclose(np.sort(a), np.sort(ra), rtol=0, atol=1e-13)
    assert np.allclose(np.sort(w), np.sort(rw), rtol=1e-11, atol=1e-12)
    # the weights are always slightly short, never long
    assert 0 < (rw.sum() - w.sum()) < 1e-11


@pytest.mark.parametrize("n", [4, 8, 15, 20])
def test_gauss_legendre_exact_for_polynomials(n):
    """Degree 2n-1 exactness on a shifted interval."""
    a, b = -0.7, 2.3
    x, w = gauss_legendre(n, a, b)
    for p in range(2 * n):
        got = float(np.sum(w * x ** p))
        want = (b ** (p + 1) - a ** (p + 1)) / (p + 1)
        # 1e-11, not machine precision: see test_gauss_legendre_matches_numpy
        assert got == pytest.approx(want, rel=1e-11, abs=1e-12)


def test_gauss_legendre_rejects_n1():
    with pytest.raises(ValueError):
        gauss_legendre(1, 0.0, 1.0)


@pytest.mark.parametrize("np_tri,degree", [(3, 2), (4, 3), (6, 4), (7, 5),
                                           (12, 6), (16, 8)])
def test_tri_quadrature_exactness(np_tri, degree):
    """Each rule integrates every monomial up to its design degree exactly on
    the reference triangle, where int xi^a eta^b = a! b! / (a+b+2)!."""
    from math import factorial

    x, y, w = tri_quadrature(np_tri)
    assert np.sum(w) == pytest.approx(1.0, rel=1e-13)
    assert np.all(x >= -1e-15) and np.all(y >= -1e-15)
    assert np.all(x + y <= 1.0 + 1e-14)
    for a in range(degree + 1):
        for b in range(degree + 1 - a):
            got = 0.5 * float(np.sum(w * x ** a * y ** b))
            want = factorial(a) * factorial(b) / factorial(a + b + 2)
            assert got == pytest.approx(want, rel=1e-12, abs=1e-14)


def test_domega_sums_to_four_pi():
    """sum(DOMEGA) = sum_j1 sin(the)*wthe * sum_j2 wphi = 2 * 2pi.

    The residual is Gauss-Legendre truncation error on int_0^pi sin, not
    round-off: it falls with NPOLE and is ~2e-14 relative at NPOLE=20.
    """
    for npole, nazim in [(10, 20), (20, 40), (32, 64)]:
        v = build_velocity_mesh(npole, nazim, 1.0)
        assert v.domega.sum() == pytest.approx(4.0 * PI, rel=1e-13)


def test_angular_moments():
    """<1> = 4pi, <cx> = <cy> = 0, <cx^2> = <cy^2>/... = 4pi Vg^2/3.

    The last identity is what makes q = -(Cv Vg^2 tau/3) grad T come out with
    the right coefficient, so it is worth asserting directly.
    """
    vg = 1.7
    v = build_velocity_mesh(24, 48, vg)
    dw = v.domega
    assert float((v.cx * dw).sum()) == pytest.approx(0.0, abs=1e-12)
    assert float((v.cy * dw).sum()) == pytest.approx(0.0, abs=1e-12)
    assert float((v.cx ** 2 * dw).sum()) == pytest.approx(4 * PI * vg ** 2 / 3, rel=1e-12)
    assert float((v.cx * v.cy * dw).sum()) == pytest.approx(0.0, abs=1e-12)


def test_nazim_forced_even():
    v = build_velocity_mesh(8, 11, 1.0)
    assert v.nazim == 10


def test_phi_covers_full_circle():
    v = build_velocity_mesh(6, 12, 1.0)
    assert v.wphi.sum() == pytest.approx(2 * PI, rel=1e-13)
    assert v.phi.min() > 0.0 and v.phi.max() < 2 * PI
