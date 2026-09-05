"""Analytic diffusion-limit reference -- from ``Out_Put_Result.f90``.

For the shipped cavity (unit square, ``T = 1`` on the north wall and ``0`` on
the other three) the Fourier limit of the Callaway model is the Laplace
solution

    T(x, y) = (2/pi) * sum_{m>=1} [((-1)^(m+1) + 1)/m] sin(m pi x)
                                   sinh(m pi y) / sinh(m pi)

truncated at 200 terms.  The even terms vanish, so this is the usual
odd-harmonic series; the mean over the square is 1/4.

The heat flux follows from Fourier's law with the Callaway conductivity
``kappa = Cv * Vg^2 * tau_R / 3`` (which is exactly the coefficient the GSIS
momentum equation carries: ``q = -(Cv*tau_R/3) grad T`` at ``Vg = 1``)::

    q = -kappa grad T

.. warning::
   ``Out_Put_Result.f90`` drops a factor of 2 when it differentiates the
   series -- it applies ``-Cv/3*TAU_R`` to the raw sum without the ``2/pi``
   prefactor's surviving ``2``.  Its ``Conduction_A.dat`` heat fluxes are
   therefore half the correct value.  ``fortran_compat=True`` reproduces that
   for byte-comparison; the default is the correct flux.  See
   docs/FORTRAN_ISSUES.md #3.
"""
from __future__ import annotations

import numpy as np

from .constants import PI

__all__ = ["fourier_temperature", "fourier_flux", "analytic_fields", "l2_error"]

N_TERMS_DEFAULT = 200


def _series(x, y, n_terms=N_TERMS_DEFAULT):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    s = np.zeros(np.broadcast(x, y).shape)
    sx = np.zeros_like(s)
    sy = np.zeros_like(s)
    if n_terms > 220:
        # sinh(m*pi) overflows binary64 beyond m ~ 226; the terms are
        # negligible long before that, so refuse rather than return inf.
        raise ValueError("n_terms > 220 overflows sinh(m*pi) in double precision")
    for m in range(1, n_terms + 1):
        a = ((-1.0) ** float(m + 1) + 1.0)
        if a == 0.0:
            continue
        mp = float(m) * PI
        den = np.sinh(mp)
        s = s + a / float(m) * np.sin(mp * x) * np.sinh(mp * y) / den
        sx = sx + a * np.cos(mp * x) * np.sinh(mp * y) / den
        sy = sy + a * np.sin(mp * x) * np.cosh(mp * y) / den
    return s, sx, sy


def fourier_temperature(x, y, n_terms=N_TERMS_DEFAULT) -> np.ndarray:
    s, _, _ = _series(x, y, n_terms)
    return s * 2.0 / PI


def fourier_flux(x, y, cv=1.0, vg=1.0, tau_r=1.0, n_terms=N_TERMS_DEFAULT,
                 fortran_compat: bool = False):
    """``(qx, qy)`` from Fourier's law with ``kappa = Cv Vg^2 tau_R / 3``."""
    _, sx, sy = _series(x, y, n_terms)
    kappa = cv * vg * vg * tau_r / 3.0
    factor = 1.0 if fortran_compat else 2.0
    return -kappa * factor * sx, -kappa * factor * sy


def analytic_fields(x, y, cv=1.0, vg=1.0, tau_r=1.0,
                    n_terms=N_TERMS_DEFAULT, fortran_compat: bool = False):
    t = fourier_temperature(x, y, n_terms)
    qx, qy = fourier_flux(x, y, cv, vg, tau_r, n_terms, fortran_compat)
    return t, qx, qy


def l2_error(sample_x, sample_y, values, weights=None, cv=1.0, vg=1.0,
             tau_r=1.0, n_terms=N_TERMS_DEFAULT) -> dict:
    """Relative L2 error of a sampled temperature field against the series.

    ``weights`` defaults to uniform, which is what a plain grid comparison
    wants; pass element areas for a mesh-based comparison.
    """
    ref = fourier_temperature(sample_x, sample_y, n_terms)
    d = np.asarray(values, dtype=float) - ref
    if weights is None:
        num = float(np.sqrt(np.mean(d * d)))
        den = float(np.sqrt(np.mean(ref * ref)))
    else:
        w = np.asarray(weights, dtype=float)
        num = float(np.sqrt(np.sum(w * d * d) / np.sum(w)))
        den = float(np.sqrt(np.sum(w * ref * ref) / np.sum(w)))
    return {"l2_abs": num, "l2_rel": num / den if den else float("nan"),
            "linf_abs": float(np.max(np.abs(d)))}
