"""Kinetic correction -- ``Correct_VDF_Calculate_Macro_Properties``.

The macroscopic solution is pushed back into the distribution through the
equilibrium and drift channels::

    tau_loc = TAU_R / Hmin(I)             (variant B: tau_loc = TAU_R)
    beta    = min(tau_loc, TAU_THR) / tau_loc

    l_T  = ( T_ACC -  T_VDF) * beta
    l_qx = (qx_ACC - qx_VDF) * beta
    l_qy = (qy_ACC - qy_VDF) * beta

    f += l_T*Cv/(4 pi) + (cx*l_qx + cy*l_qy) * (TAU_C/TAU_N) * 3/(4 pi Vg^2)

``beta`` damps the correction where the local element is optically thin, so
the acceleration cannot overshoot in the ballistic layer.  Note the heat-flux
channel enters through ``TAU_C/TAU_N``: with the shipped ``TAU_N = 1e5`` it is
numerically negligible, and the correction is effectively a temperature
update -- but ``T_s``/``Qx_s``/``Qy_s`` are still set to the *full* corrected
macroscopic values, which is what feeds the next sweep.
"""
from __future__ import annotations

from .._numba import njit, prange

__all__ = ["correct_vdf"]


@njit(cache=True, parallel=True)
def _correct_kernel(vdf, uq, cxv, cyv, domega, int_tri, tri_hmin,
                    cv, vg, tau_r, tau_n, tau_c, tau_thr, pi, nd,
                    use_hmin, ts, qxs, qys, temp, qx, qy, defect, mkin):
    n_tris = vdf.shape[1]
    ndir = vdf.shape[0]
    drift = tau_c / tau_n * 3.0 / 4.0 / pi / vg / vg
    for i in prange(n_tris):
        if use_hmin:
            tau_loc = tau_r / tri_hmin[i]
        else:
            tau_loc = tau_r
        beta = min(tau_loc, tau_thr) / tau_loc
        temp[i] = 0.0
        qx[i] = 0.0
        qy[i] = 0.0
        for m in range(nd):
            t_vdf = 0.0
            qx_vdf = 0.0
            qy_vdf = 0.0
            for d in range(ndir):
                f = vdf[d, i, m]
                t_vdf += f * domega[d]
                qx_vdf += cxv[d] * f * domega[d]
                qy_vdf += cyv[d] * f * domega[d]
            t_vdf = t_vdf / cv
            mkin[m, i] = t_vdf
            mkin[nd + m, i] = qx_vdf
            mkin[2 * nd + m, i] = qy_vdf

            # ``defect`` is zero for the shipped scheme, in which case these
            # three lines are the damped blend unchanged.
            l_t = (uq[m, i] - defect[m, i] - t_vdf) * beta
            l_qx = (uq[nd + m, i] - defect[nd + m, i] - qx_vdf) * beta
            l_qy = (uq[2 * nd + m, i] - defect[2 * nd + m, i] - qy_vdf) * beta

            eq = l_t * cv / 4.0 / pi
            for d in range(ndir):
                vdf[d, i, m] = vdf[d, i, m] + eq + (cxv[d] * l_qx + cyv[d] * l_qy) * drift

            ts[m, i] = l_t + t_vdf
            qxs[m, i] = l_qx + qx_vdf
            qys[m, i] = l_qy + qy_vdf

            temp[i] += (l_t + t_vdf) * int_tri[m, i]
            qx[i] += (l_qx + qx_vdf) * int_tri[m, i]
            qy[i] += (l_qy + qy_vdf) * int_tri[m, i]


def correct_vdf(vdf, uq, cxv, cyv, domega, int_tri, tri_hmin, cv, vg,
                tau_r, tau_n, tau_c, tau_thr, pi, ndof_tri, mom,
                use_hmin: bool = True, defect=None, mkin=None):
    """``defect``/``mkin`` are ``(3*ndof_tri, n_tris)``.

    ``defect`` is subtracted from the macroscopic solution before the blend;
    passing zeros (or ``None``) gives the shipped blend exactly.
    ``mkin`` receives the kinetic moments ``T, qx, qy`` per DOF, which the
    caller needs to update the defect and which are computed here anyway.
    """
    import numpy as _np

    if defect is None:
        defect = _np.zeros((3 * ndof_tri, vdf.shape[1]))
    if mkin is None:
        mkin = _np.empty((3 * ndof_tri, vdf.shape[1]))
    _correct_kernel(vdf, uq, cxv, cyv, domega, int_tri, tri_hmin,
                    cv, vg, tau_r, tau_n, tau_c, tau_thr, pi, ndof_tri,
                    use_hmin, mom.ts, mom.qxs, mom.qys,
                    mom.temp, mom.qx, mom.qy, defect, mkin)
    return mom
