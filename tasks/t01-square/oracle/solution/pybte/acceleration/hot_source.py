"""Higher-order-term source -- ``Calculate_SRC_ACC_HoTfromDVM``.

The macroscopic system is *not* a closure approximation: the stress terms it
carries are computed from the kinetic solution, which is what makes GSIS
converge to the same discrete fixed point as CIS rather than to a moment
model.  This module builds that coupling.

For each element and each element DOF ``M`` the kinetic distribution is
contracted with the x- and y-derivative tensors and with the third-order
angular moments

    PIxx  <-  cx(5cx^2-3) d_x f  +  cy(5cx^2-1) d_y f
    PIxy  <-  cy(5cx^2-1) d_x f  +  cx(5cy^2-1) d_y f
    PIyy  <-  cx(5cy^2-1) d_x f  +  cy(5cy^2-3) d_y f

and the four stress rows of ``AA_SRC`` follow.  The whole vector is then
premultiplied by ``inv(AA_SOL)``, so ``AA_SRC`` as stored is already the
particular solution of the local problem.
"""
from __future__ import annotations

import numpy as np

from .._numba import njit, prange

__all__ = ["hot_source"]


@njit(cache=True, parallel=True)
def _hot_kernel(vdf, gx, gy, cxv, cyv, domega, tc5, nd, out):
    """``out`` is ``(n_tris, 7*nd)``; only the four stress blocks are filled."""
    n_tris = vdf.shape[1]
    ndir = vdf.shape[0]
    for i in prange(n_tris):
        for m in range(nd):
            pixx = 0.0
            pixy = 0.0
            piyy = 0.0
            for d in range(ndir):
                cx = cxv[d]
                cy = cyv[d]
                dw = domega[d]
                t = 0.0
                for l in range(nd):
                    t += vdf[d, i, l] * gx[i, l, m]
                pixx += cx * (5.0 * cx * cx - 3.0) * t * dw
                pixy += cy * (5.0 * cx * cx - 1.0) * t * dw
                piyy += cx * (5.0 * cy * cy - 1.0) * t * dw
                t = 0.0
                for l in range(nd):
                    t += vdf[d, i, l] * gy[i, l, m]
                pixx += cy * (5.0 * cx * cx - 1.0) * t * dw
                pixy += cx * (5.0 * cy * cy - 1.0) * t * dw
                piyy += cy * (5.0 * cy * cy - 3.0) * t * dw
            out[i, 3 * nd + m] = (pixx + 0.5 * piyy) * tc5
            out[i, 4 * nd + m] = 0.5 * pixy * tc5
            out[i, 5 * nd + m] = 0.5 * pixy * tc5
            out[i, 6 * nd + m] = (0.5 * pixx + piyy) * tc5


def hot_source(vdf, gx, gy, cxv, cyv, domega, tau_c, ndof_tri, inv_aa_sol,
               scale=1.0, work=None):
    """Return ``AA_SRC`` with shape ``(7*ndof_tri, n_tris)`` (Fortran layout).

    ``scale`` is ``TAU_R`` for acceleration variant B and 1 for variant A.
    """
    n_tris = vdf.shape[1]
    if work is None or work.shape != (n_tris, 7 * ndof_tri):
        work = np.zeros((n_tris, 7 * ndof_tri))
    else:
        work[:] = 0.0
    _hot_kernel(vdf, gx, gy, cxv, cyv, domega,
                tau_c * scale / 5.0, ndof_tri, work)
    out = np.einsum("ipq,iq->pi", inv_aa_sol, work, optimize=True)
    return np.ascontiguousarray(out), work
