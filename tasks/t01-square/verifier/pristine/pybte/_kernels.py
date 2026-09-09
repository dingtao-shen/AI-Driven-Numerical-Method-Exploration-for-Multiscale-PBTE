"""Numba kernels: dense LU, the transport sweep, and the moment reduction.

Two things in here need explaining.

**Why hand-rolled LU.**  The Fortran calls ``DGETRF``/``DGETRS`` per element
per direction per iteration.  Those factorisations depend only on geometry,
direction and ``TAU_C``, all iteration-invariant, so we hoist them out of the
loop.  A batched LAPACK call would be ideal but scipy has none,
and 160k tiny factorisations through the Python layer costs more than the
sweep itself.  The kernels below reproduce LAPACK's unblocked ``DGETF2`` and
``DTRSM`` *operation for operation*:

* pivot = first index of maximum |a| in the column (``IDAMAX``);
* the multipliers are formed as ``a * (1/pivot)``, not ``a / pivot``, which is
  what ``DGETF2`` does whenever ``|pivot| >= SFMIN``;
* the rank-1 update is a single fused ``a[i][j] -= a[i][k]*a[k][j]`` per entry,
  so there is no summation order to get wrong;
* back substitution accumulates from ``k = n-1`` **downwards**, matching
  ``DTRSM``'s column sweep; the naive increasing-``k`` dot product would give
  a different last bit.

The result is bit-identical to the Fortran's LAPACK for every case tested.

**Why the direction loop is the parallel one.**  Directions are completely
independent -- each writes its own slice of ``vdf`` and reads only its own --
so ``prange`` over directions is safe and mirrors the Fortran's OpenMP.  The
sweep *within* a direction is strictly sequential: it reads neighbour values
that earlier elements in the same pass have already overwritten.  That
Gauss-Seidel coupling is the whole point of the ordering and must not be
parallelised away.
"""
from __future__ import annotations

import numpy as np

from ._numba import njit, prange

FOUR_PI = 4.0


# ---------------------------------------------------------------------------
# LAPACK look-alikes
# ---------------------------------------------------------------------------
@njit(cache=True, inline="always")
def lu_factor_inplace(a, ipiv):
    """``DGETF2`` on a row-major ``(n, n)`` block.  Returns 0 or the failing
    column (1-based), like ``INFO``."""
    n = a.shape[0]
    info = 0
    for k in range(n):
        p = k
        big = abs(a[k, k])
        for i in range(k + 1, n):
            v = abs(a[i, k])
            if v > big:
                big = v
                p = i
        ipiv[k] = p
        if a[p, k] == 0.0:
            if info == 0:
                info = k + 1
            continue
        if p != k:
            for j in range(n):
                t = a[k, j]
                a[k, j] = a[p, j]
                a[p, j] = t
        r = 1.0 / a[k, k]
        for i in range(k + 1, n):
            a[i, k] = a[i, k] * r
        for i in range(k + 1, n):
            aik = a[i, k]
            if aik != 0.0:
                for j in range(k + 1, n):
                    a[i, j] -= aik * a[k, j]
    return info


@njit(cache=True, inline="always")
def lu_solve_inplace(a, ipiv, b):
    """``DGETRS('N')`` for a single right-hand side."""
    n = a.shape[0]
    for k in range(n):
        p = ipiv[k]
        if p != k:
            t = b[k]
            b[k] = b[p]
            b[p] = t
    # L y = Pb   (unit lower)
    for k in range(n):
        bk = b[k]
        if bk != 0.0:
            for i in range(k + 1, n):
                b[i] -= bk * a[i, k]
    # U x = y    -- accumulate downwards, as DTRSM does
    for k in range(n - 1, -1, -1):
        b[k] = b[k] / a[k, k]
        bk = b[k]
        if bk != 0.0:
            for i in range(k):
                b[i] -= bk * a[i, k]


# ---------------------------------------------------------------------------
# element operator A_SOL
# ---------------------------------------------------------------------------
@njit(cache=True, inline="always")
def _fill_a_sol(a, i, cxd, cyd, inv_tau_c, mass, gx, gy, fcm, nx, ny):
    """``A_SOL(M,L)`` for element ``i`` and one direction, exactly as
    ``Solvers.f90`` writes it."""
    nd = a.shape[0]
    s0 = cxd * nx[i, 0] + cyd * ny[i, 0]
    s1 = cxd * nx[i, 1] + cyd * ny[i, 1]
    s2 = cxd * nx[i, 2] + cyd * ny[i, 2]
    w0 = 0.5 * (s0 + abs(s0))
    w1 = 0.5 * (s1 + abs(s1))
    w2 = 0.5 * (s2 + abs(s2))
    for l in range(nd):
        for m in range(nd):
            v = (inv_tau_c * mass[m, l, i]
                 - cxd * gx[i, m, l] - cyd * gy[i, m, l])
            v = v + w0 * fcm[i, 0, m, l]
            v = v + w1 * fcm[i, 1, m, l]
            v = v + w2 * fcm[i, 2, m, l]
            a[m, l] = v


@njit(cache=True, parallel=True)
def build_operators(order, cxv, cyv, inv_tau_c, mass, gx, gy, fcm, nx, ny,
                    lu, piv):
    """Factorise ``A_SOL`` once for every ``(direction, element)`` pair."""
    ndir = cxv.shape[0]
    n_tris = mass.shape[2]
    nd = mass.shape[0]
    bad = 0
    for d in prange(ndir):
        a = np.empty((nd, nd))
        p = np.empty(nd, dtype=np.int32)
        for i in range(n_tris):
            _fill_a_sol(a, i, cxv[d], cyv[d], inv_tau_c, mass, gx, gy, fcm, nx, ny)
            info = lu_factor_inplace(a, p)
            if info != 0:
                bad += 1
            for m in range(nd):
                piv[d, i, m] = p[m]
                for l in range(nd):
                    lu[d, i, m, l] = a[m, l]
    return bad


# ---------------------------------------------------------------------------
# the sweep
# ---------------------------------------------------------------------------
@njit(cache=True, inline="always")
def _fill_a_src(src, i, cxd, cyd, cv, vg, tau_r, tau_n, pi,
                mass, ts, qxs, qys):
    """The scattering source, in ``Solvers.f90``'s exact expression order."""
    nd = src.shape[0]
    for m in range(nd):
        acc = 0.0
        for l in range(nd):
            mlm = mass[l, m, i]
            acc = acc + cv * ts[l, i] / 4.0 / pi / tau_r * mlm
            acc = acc + 1.0 / tau_n * mlm * (
                cv * ts[l, i] / 4.0 / pi
                + 3.0 / 4.0 / pi * (qxs[l, i] * cxd + qys[l, i] * cyd) / vg / vg)
        src[m] = acc


@njit(cache=True, inline="always")
def _add_face_sources(src, i, d, cxd, cyd, cv, pi,
                      ttfc, nx, ny, tri_faces, face_bc, neighbour,
                      periodic_tri, bc_type, bc_temp, flux_wall, vdf):
    """Upwind inflow, dispatching on the boundary type."""
    nd = src.shape[0]
    for il in range(3):
        speed = cxd * nx[i, il] + cyd * ny[i, il]
        w = 0.5 * (speed - abs(speed))
        if w == 0.0:
            continue
        fc = tri_faces[i, il]
        bc = face_bc[fc]
        if bc < 0:
            iext = neighbour[i, il]
            for l in range(nd):
                fl = vdf[d, iext, l]
                for m in range(nd):
                    src[m] -= (w * ttfc[i, il, m, l]) * fl
        else:
            typ = bc_type[bc]
            if typ == 3:
                iext = periodic_tri[fc]
                for l in range(nd):
                    fl = vdf[d, iext, l]
                    for m in range(nd):
                        src[m] -= (w * ttfc[i, il, m, l]) * fl
            elif typ == 2:
                for m in range(nd):
                    src[m] -= w * (-flux_wall[fc, m])
            else:
                tw = cv / 4.0 / pi * bc_temp[bc]
                for m in range(nd):
                    src[m] -= w * (tw * ttfc[i, il, m, 0])


@njit(cache=True, parallel=True)
def sweep_precomputed(order, cxv, cyv, cv, vg, tau_r, tau_n, pi,
                      mass, ttfc, nx, ny, tri_faces, face_bc, neighbour,
                      periodic_tri, bc_type, bc_temp, flux_wall,
                      ts, qxs, qys, lu, piv, vdf):
    """One transport sweep using pre-factorised element operators."""
    ndir = cxv.shape[0]
    n_tris = order.shape[1]
    nd = mass.shape[0]
    for d in prange(ndir):
        src = np.empty(nd)
        cxd = cxv[d]
        cyd = cyv[d]
        for k in range(n_tris):
            i = order[d, k]
            _fill_a_src(src, i, cxd, cyd, cv, vg, tau_r, tau_n, pi, mass, ts, qxs, qys)
            _add_face_sources(src, i, d, cxd, cyd, cv, pi, ttfc, nx, ny,
                              tri_faces, face_bc, neighbour, periodic_tri,
                              bc_type, bc_temp, flux_wall, vdf)
            lu_solve_inplace(lu[d, i], piv[d, i], src)
            for m in range(nd):
                vdf[d, i, m] = src[m]


@njit(cache=True, parallel=True)
def sweep_onthefly(order, cxv, cyv, cv, vg, tau_r, tau_n, tau_c, pi,
                   mass, gx, gy, fcm, ttfc, nx, ny, tri_faces, face_bc,
                   neighbour, periodic_tri, bc_type, bc_temp, flux_wall,
                   ts, qxs, qys, vdf):
    """One transport sweep, re-factorising ``A_SOL`` exactly as the Fortran
    does.  Slower, but has no ``O(N_TRIS*NDIR*NDOF^2)`` memory cost."""
    ndir = cxv.shape[0]
    n_tris = order.shape[1]
    nd = mass.shape[0]
    inv_tau_c = 1.0 / tau_c
    for d in prange(ndir):
        src = np.empty(nd)
        a = np.empty((nd, nd))
        p = np.empty(nd, dtype=np.int32)
        cxd = cxv[d]
        cyd = cyv[d]
        for k in range(n_tris):
            i = order[d, k]
            _fill_a_sol(a, i, cxd, cyd, inv_tau_c, mass, gx, gy, fcm, nx, ny)
            _fill_a_src(src, i, cxd, cyd, cv, vg, tau_r, tau_n, pi, mass, ts, qxs, qys)
            _add_face_sources(src, i, d, cxd, cyd, cv, pi, ttfc, nx, ny,
                              tri_faces, face_bc, neighbour, periodic_tri,
                              bc_type, bc_temp, flux_wall, vdf)
            lu_factor_inplace(a, p)
            lu_solve_inplace(a, p, src)
            for m in range(nd):
                vdf[d, i, m] = src[m]


# ---------------------------------------------------------------------------
# residual of the discrete transport system
# ---------------------------------------------------------------------------
@njit(cache=True, parallel=True)
def transport_residual(cxv, cyv, domega, cv, vg, tau_r, tau_n, tau_c, pi,
                       mass, gx, gy, fcm, ttfc, nx, ny, tri_faces, face_bc,
                       neighbour, periodic_tri, bc_type, bc_temp, flux_wall,
                       ts, qxs, qys, vdf, out_num, out_den):
    """``||A_SOL f - A_SRC||`` with *no* sweep: a genuine residual of the
    coupled discrete system, zero exactly at the fixed point."""
    ndir = cxv.shape[0]
    n_tris = mass.shape[2]
    nd = mass.shape[0]
    inv_tau_c = 1.0 / tau_c
    for d in prange(ndir):
        src = np.empty(nd)
        a = np.empty((nd, nd))
        num = 0.0
        den = 0.0
        cxd = cxv[d]
        cyd = cyv[d]
        for i in range(n_tris):
            _fill_a_sol(a, i, cxd, cyd, inv_tau_c, mass, gx, gy, fcm, nx, ny)
            _fill_a_src(src, i, cxd, cyd, cv, vg, tau_r, tau_n, pi, mass, ts, qxs, qys)
            _add_face_sources(src, i, d, cxd, cyd, cv, pi, ttfc, nx, ny,
                              tri_faces, face_bc, neighbour, periodic_tri,
                              bc_type, bc_temp, flux_wall, vdf)
            for m in range(nd):
                r = -src[m]
                for l in range(nd):
                    r += a[m, l] * vdf[d, i, l]
                num += r * r * domega[d]
                den += src[m] * src[m] * domega[d]
        out_num[d] = num
        out_den[d] = den


# ---------------------------------------------------------------------------
# moments
# ---------------------------------------------------------------------------
@njit(cache=True)
def moments_kernel(vdf, cxv, cyv, domega, int_tri, cv, ts, qxs, qys,
                   temp, qx, qy):
    """``Calculate_Macro_Properties``."""
    ndir = vdf.shape[0]
    n_tris = vdf.shape[1]
    nd = vdf.shape[2]
    for i in range(n_tris):
        temp[i] = 0.0
        qx[i] = 0.0
        qy[i] = 0.0
    for i in range(n_tris):
        for l in range(nd):
            ss = 0.0
            ssx = 0.0
            ssy = 0.0
            for d in range(ndir):
                f = vdf[d, i, l]
                ss += f * domega[d]
                ssx += cxv[d] * f * domega[d]
                ssy += cyv[d] * f * domega[d]
            ts[l, i] = ss / cv
            qxs[l, i] = ssx
            qys[l, i] = ssy
            temp[i] += ss * int_tri[l, i] / cv
            qx[i] += ssx * int_tri[l, i]
            qy[i] += ssy * int_tri[l, i]


@njit(cache=True)
def sum_kernel(a):
    """Sequential sum, matching Fortran's ``SUM`` (numpy's is pairwise)."""
    s = 0.0
    for i in range(a.shape[0]):
        s += a[i]
    return s


@njit(cache=True)
def residual_kernel(temp, t_old):
    """``Calculate_Residual_T``, accumulated element by element in index
    order so that it matches the Fortran bit for bit (numpy's pairwise
    summation would not)."""
    ss1 = 0.0
    ss2 = 0.0
    for i in range(temp.shape[0]):
        d = temp[i] - t_old[i]
        ss1 += d * d
        ss2 += temp[i] * temp[i]
    return np.sqrt(ss1 / ss2)


@njit(cache=True, parallel=True)
def wall_flux_kernel(vdf, cxv, cyv, domega, fcm, nx, ny,
                     face_bc, face_tri, face_lfc, bc_type, flux_wall):
    """``Calculate_FLUX_WALL`` -- the diffuse (adiabatic) wall emission that
    makes the net normal energy flux vanish.  Dead code in the Fortran; wired
    in here."""
    n_faces = face_bc.shape[0]
    ndir = cxv.shape[0]
    nd = fcm.shape[2]
    for f in prange(n_faces):
        bc = face_bc[f]
        if bc < 0 or bc_type[bc] != 2:
            continue
        if face_tri[f, 0] < 0:
            tri = face_tri[f, 1]
            lfc = face_lfc[f, 1]
        else:
            tri = face_tri[f, 0]
            lfc = face_lfc[f, 0]
        n1 = -nx[tri, lfc]
        n2 = -ny[tri, lfc]
        for m in range(nd):
            ss1 = 0.0
            ss2 = 0.0
            for d in range(ndir):
                u = cxv[d] * n1 + cyv[d] * n2
                if u < 0.0:
                    acc = 0.0
                    for l in range(nd):
                        acc += vdf[d, tri, l] * fcm[tri, lfc, m, l]
                    ss1 += u * domega[d] * acc
                else:
                    ss2 -= u * domega[d]
            flux_wall[f, m] = -ss1 / ss2
