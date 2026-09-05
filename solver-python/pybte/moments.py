"""Macroscopic moments and residuals -- from ``Velocity_Distribution.f90``.

``Calculate_Macro_Properties``::

    Cv * T_s(L,I) = sum_dir  VDF(L,I,dir) * DOMEGA
    qx_s(L,I)     = sum_dir  CX * VDF * DOMEGA
    qy_s(L,I)     = sum_dir  CY * VDF * DOMEGA
    Temp(I)       = sum_L  Cv*T_s(L,I) * INT_NODFUNC_TRI(L,I) / Cv

``Temp``/``Qx``/``Qy`` are element *integrals*, not averages: summed over the
mesh they give the domain integral of the field (0.25 for the shipped cavity,
which is the mean of the Laplace solution over the unit square).

``Calculate_Residual_T`` is the Fortran's stopping criterion and is kept
verbatim as ``residual_iterate``; see §7.6 for the true residual that is
reported alongside it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._kernels import moments_kernel, residual_kernel, sum_kernel

__all__ = ["Moments", "compute_moments", "residual_iterate"]


@dataclass
class Moments:
    """Per-DOF and per-element macroscopic fields."""
    ts: np.ndarray     # (ndof_tri, n_tris)  perturbed temperature
    qxs: np.ndarray    # (ndof_tri, n_tris)
    qys: np.ndarray    # (ndof_tri, n_tris)
    temp: np.ndarray   # (n_tris,)  element integral of T
    qx: np.ndarray     # (n_tris,)
    qy: np.ndarray     # (n_tris,)

    @property
    def mass(self) -> float:
        return float(sum_kernel(self.temp))

    @classmethod
    def zeros(cls, ndof_tri: int, n_tris: int) -> "Moments":
        return cls(ts=np.zeros((ndof_tri, n_tris)),
                   qxs=np.zeros((ndof_tri, n_tris)),
                   qys=np.zeros((ndof_tri, n_tris)),
                   temp=np.zeros(n_tris), qx=np.zeros(n_tris), qy=np.zeros(n_tris))


def compute_moments(vdf, cxv, cyv, domega, int_tri, cv, mom: Moments,
                    zero_qy: bool = True) -> Moments:
    """In-place update of ``mom`` from ``vdf`` (shape ``(ndir, n_tris, ndof)``)."""
    moments_kernel(vdf, cxv, cyv, domega, int_tri, cv,
                   mom.ts, mom.qxs, mom.qys, mom.temp, mom.qx, mom.qy, zero_qy)
    return mom


def residual_iterate(temp: np.ndarray, t_old: np.ndarray) -> float:
    """``Calculate_Residual_T``: ``sqrt(sum (T-T_old)^2 / sum T^2)``.

    Note what this is *not*: it is the normalised change between successive
    iterates -- a *step*, not an *error*.  For a linearly converging sequence
    the remaining error is the sum of all future steps,
    ``r * rho/(1 - rho)``, and ``rho -> 1`` as the medium becomes optically
    thick.  On the shipped mesh at ``tau_R = 1e-4`` CIS reports ``2.5e-6``
    after 200 000 iterations while its answer is 89% wrong.

    That pseudo-convergence is exactly the phenomenon the benchmark studies,
    so this criterion is preserved verbatim as the stopping rule.
    ``RunRecord.error_estimate`` reports the implied error alongside it.
    """
    return float(residual_kernel(temp, t_old))
