"""Reference transport sweep: plain numpy + LAPACK, no numba, no precompute.

This is the honest transliteration of ``Solvers.f90`` -- it rebuilds and
re-factorises ``A_SOL`` for every element, every direction, every call, and it
calls the real ``DGETRF``/``DGETRS`` through :mod:`scipy.linalg`.  It is far
too slow for production (minutes per sweep on the shipped case) and exists for
exactly one reason: to pin down what :mod:`pybte._kernels` must reproduce.

``tests/unit/test_sweep_kernels.py`` asserts that this and the jitted sweep
agree, so an optimisation that quietly changes the arithmetic cannot pass.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import lu_factor, lu_solve

from .constants import PI

__all__ = ["sweep_reference"]


def sweep_reference(ctx, mom, vdf) -> None:
    """One sweep, in place on ``vdf`` with shape ``(ndir, n_tris, ndof)``."""
    nd = ctx.mass.shape[0]
    n_tris = ctx.order.shape[1]
    ndir = ctx.cxv.size
    cv, vg = ctx.cv, ctx.vg
    tau_r, tau_n, tau_c = ctx.tau_r, ctx.tau_n, ctx.tau_c

    for d in range(ndir):
        cxd, cyd = ctx.cxv[d], ctx.cyv[d]
        for k in range(n_tris):
            i = int(ctx.order[d, k])

            # ---- A_SOL ------------------------------------------------
            a = (1.0 / tau_c) * ctx.mass[:, :, i] \
                - cxd * ctx.gx[i] - cyd * ctx.gy[i]
            for il in range(3):
                speed = cxd * ctx.nx[i, il] + cyd * ctx.ny[i, il]
                a = a + 0.5 * (speed + abs(speed)) * ctx.fcm[i, il]

            # ---- A_SRC ------------------------------------------------
            # Accumulated in the Fortran's exact order: both terms per L
            # inside one loop over L, not two separate matvecs.  The two
            # differ in the last bit and the difference compounds over
            # millions of iterations.
            ts = mom.ts[:, i]
            qxs = mom.qxs[:, i]
            qys = mom.qys[:, i]
            src = np.zeros(nd)
            for m in range(nd):
                acc = 0.0
                for l in range(nd):
                    mlm = ctx.mass[l, m, i]
                    acc = acc + cv * ts[l] / 4.0 / PI / tau_r * mlm
                    acc = acc + 1.0 / tau_n * mlm * (
                        cv * ts[l] / 4.0 / PI
                        + 3.0 / 4.0 / PI * (qxs[l] * cxd + qys[l] * cyd) / vg / vg)
                src[m] = acc

            for il in range(3):
                speed = cxd * ctx.nx[i, il] + cyd * ctx.ny[i, il]
                w = 0.5 * (speed - abs(speed))
                if w == 0.0:
                    continue
                fc = int(ctx.tri_faces[i, il])
                bc = int(ctx.face_bc[fc])
                if bc < 0 or ctx.bc_type[bc] == 3:
                    iext = (int(ctx.neighbour[i, il]) if bc < 0
                            else int(ctx.periodic_tri[fc]))
                    for l in range(nd):
                        src -= (w * ctx.ttfc[i, il, :, l]) * vdf[d, iext, l]
                elif ctx.bc_type[bc] == 2:
                    src -= w * (-ctx.flux_wall[fc])
                else:
                    tw = cv / 4.0 / PI * ctx.bc_temp[bc]
                    src -= w * (tw * ctx.ttfc[i, il, :, 0])

            lu, piv = lu_factor(a)
            vdf[d, i] = lu_solve((lu, piv), src)
