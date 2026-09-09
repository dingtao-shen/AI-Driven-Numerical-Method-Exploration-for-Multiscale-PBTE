"""Quadrature rules -- the port of ``USD_Math.f90``.

Two things live here:

``gauss_legendre``
    A *literal* transcription of the Fortran's hand-rolled Newton iteration,
    not a call to :func:`numpy.polynomial.legendre.leggauss`.  The nodes agree
    with numpy's to ~1e-16, but the Fortran's own iteration is reproduced so
    that no last-bit difference can leak into the angular quadrature and,
    through it, into the iteration counts we are trying to match exactly.

``tri_quadrature``
    The symmetric triangle rules for NP = 3, 4, 6, 7, 12, 16 on the reference
    triangle (0,0)-(1,0)-(0,1), weights normalised to sum to 1.
"""
from __future__ import annotations

import numpy as np

from .constants import PI

__all__ = ["gauss_legendre", "tri_quadrature", "leggauss_reference"]


_GL_CORE_CACHE: dict = {}


def gauss_legendre_core(n: int, tol: float = 1.0e-13, max_sweeps: int = 200):
    """The interval-independent part of ``GaussLegendre``: ``(y, lp)``.

    The Newton iteration in the Fortran never touches ``a`` or ``b``, so the
    roots and derivatives can be cached and mapped afterwards.  The mapping
    below is written exactly as the Fortran writes it, so this is a pure
    speed-up with no change to a single bit.
    """
    key = (n, tol)
    if key in _GL_CORE_CACHE:
        return _GL_CORE_CACHE[key]
    y, lp = _gl_newton(n, tol, max_sweeps)
    _GL_CORE_CACHE[key] = (y, lp)
    return y, lp


def gauss_legendre(n: int, a: float, b: float, tol: float = 1.0e-13,
                   max_sweeps: int = 200):
    """Gauss-Legendre abscissae and weights on ``[a, b]``.

    Transcribes ``USD_Math.f90::GaussLegendre`` statement for statement,
    including the ``0.27/N1*sin(...)`` perturbation of the initial guess and
    the ``N2/N1`` scaling of the derivative that cancels in the weights.

    Returns
    -------
    (absc, weig) : tuple of ndarray, shape (n,)
    """
    y, lp = gauss_legendre_core(n, tol, max_sweeps)
    n1, n2 = n, n + 1
    absc = (a * (1.0 - y) + b * (1.0 + y)) / 2.0
    weig = n2 * n2 * (b - a) / ((1.0 - y * y) * lp * lp) / n1 / n1
    return absc, weig


def _gl_newton(n: int, tol: float, max_sweeps: int):
    if n < 2:
        raise ValueError("the Fortran initial guess divides by (N-1); n >= 2 required")

    n1 = n
    n2 = n + 1

    i = np.arange(1, n1 + 1, dtype=np.float64)
    y = (np.cos((2.0 * (i - 1.0) + 1.0) * PI / (2.0 * (n - 1.0) + 2.0))
         + 0.27 / n1 * np.sin(PI * (-1.0 + i * 2.0 / (n1 - 1.0)) * (n - 1.0) / n2))

    lp = np.zeros(n1)
    leg = np.zeros((n1, n2))

    for _ in range(max_sweeps):
        leg[:, 0] = 1.0
        leg[:, 1] = y
        for k in range(2, n1 + 1):
            leg[:, k] = ((2.0 * k - 1.0) * y * leg[:, k - 1] - (k - 1) * leg[:, k - 2]) / k
        lp = n2 * (leg[:, n1 - 1] - y * leg[:, n2 - 1]) / (1.0 - y * y)
        y0 = y
        y = y0 - leg[:, n2 - 1] / lp
        if np.max(np.abs(y - y0)) < tol:
            break
    else:  # pragma: no cover - the Fortran loops forever instead
        raise RuntimeError("gauss_legendre did not converge")

    return y, lp


def leggauss_reference(n: int, a: float, b: float):
    """Independent cross-check via numpy (used only in the unit tests)."""
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (a * (1.0 - x) + b * (1.0 + x)), 0.5 * (b - a) * w


# ---------------------------------------------------------------------------
# Triangle rules.  Ordering matters: ``NODFUN_QUA_P`` is indexed by the
# quadrature point, so the point order is part of the reproducible state.
# ---------------------------------------------------------------------------
def tri_quadrature(np_tri: int):
    """Return ``(abscissa_x, abscissa_y, weight)`` for the reference triangle.

    Transcribes ``USD_Math.f90::TRI_QUADRATURE``.  Weights sum to 1; the
    physical integral is ``2*area * sum(w * f)``.
    """
    x = np.zeros(np_tri)
    y = np.zeros(np_tri)
    w = np.zeros(np_tri)

    if np_tri == 3:
        w[:] = 1.0 / 3.0
        x[0], y[0] = 0.5, 0.0
        x[1], y[1] = 0.5, 0.5
        x[2], y[2] = 0.0, 0.5
    elif np_tri == 4:
        w[0] = -27.0 / 48.0
        w[1:4] = 25.0 / 48.0
        x[0], y[0] = 1.0 / 3.0, 1.0 / 3.0
        x[1], y[1] = 0.2, 0.2
        x[2], y[2] = 0.6, 0.2
        x[3], y[3] = 0.2, 0.6
    elif np_tri == 6:
        w[0:3] = 0.109951743655322
        w[3:6] = 0.223381589678011
        x[0], y[0] = 0.091576213509771, 0.091576213509771
        x[1], y[1] = 0.816847572980459, 0.091576213509771
        x[2], y[2] = 0.091576213509771, 0.816847572980459
        x[3], y[3] = 0.445948490915965, 0.445948490915965
        x[4], y[4] = 0.108103018168070, 0.445948490915965
        x[5], y[5] = 0.445948490915965, 0.108103018168070
    elif np_tri == 7:
        w[0] = 0.225
        w[1:4] = 0.125939180544827
        w[4:7] = 0.132394152788506
        x[0], y[0] = 1.0 / 3.0, 1.0 / 3.0
        x[1], y[1] = 0.101286507323456, 0.101286507323456
        x[2], y[2] = 0.797426985353087, 0.101286507323456
        x[3], y[3] = 0.101286507323456, 0.797426985353087
        x[4], y[4] = 0.470142064105115, 0.059715871789770
        x[5], y[5] = 0.470142064105115, 0.470142064105115
        x[6], y[6] = 0.059715871789770, 0.470142064105115
    elif np_tri == 12:
        w[0:3] = 0.050844906370207
        w[3:6] = 0.116786275726379
        w[6:12] = 0.082851075618374
        a = 0.063089014491502
        x[0], y[0] = a, 1.0 - 2.0 * a
        x[1], y[1] = a, a
        x[2], y[2] = 1.0 - 2.0 * a, a
        a = 0.249286745170911
        x[3], y[3] = a, 1.0 - 2.0 * a
        x[4], y[4] = a, a
        x[5], y[5] = 1.0 - 2.0 * a, a
        a = 0.636502499121399
        b = 0.310352451033785
        x[6], y[6] = b, 1.0 - a - b
        x[7], y[7] = 1.0 - a - b, b
        x[8], y[8] = a, b
        x[9], y[9] = b, a
        x[10], y[10] = 1.0 - a - b, a
        x[11], y[11] = a, 1.0 - a - b
    elif np_tri == 16:
        w[0] = 0.144315607677787
        w[1:4] = 0.103217370534718
        w[4:7] = 0.032458497623198
        w[7:10] = 0.095091634267285
        w[10:16] = 0.027230314174435
        x[0], y[0] = 1.0 / 3.0, 1.0 / 3.0
        a = 0.170569307751760
        x[1], y[1] = a, 1.0 - 2.0 * a
        x[2], y[2] = a, a
        x[3], y[3] = 1.0 - 2.0 * a, a
        a = 0.050547228317031
        x[4], y[4] = a, 1.0 - 2.0 * a
        x[5], y[5] = a, a
        x[6], y[6] = 1.0 - 2.0 * a, a
        a = 0.459292588292723
        x[7], y[7] = a, 1.0 - 2.0 * a
        x[8], y[8] = a, a
        x[9], y[9] = 1.0 - 2.0 * a, a
        a = 0.263112829634638
        b = 0.008394777409958
        x[10], y[10] = b, 1.0 - a - b
        x[11], y[11] = 1.0 - a - b, b
        x[12], y[12] = a, b
        x[13], y[13] = b, a
        x[14], y[14] = 1.0 - a - b, a
        x[15], y[15] = a, 1.0 - a - b
    else:
        raise ValueError(f"no triangle rule for NP_TRI={np_tri}")

    return x, y, w
