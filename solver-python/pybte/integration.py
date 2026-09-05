"""Precomputed integral tensors -- the port of ``Integration.f90``.

Every integral the solver ever needs is assembled once here.  Two techniques
are used, and which one is used for which tensor is part of the reference
behaviour we have to reproduce:

* **analytic** -- contract the nodal-basis monomial coefficients with the
  closed-form monomial moments ``int_T xi^a eta^b = a! b! / (a+b+2)!``
  (and the edge analogues).  Used for the mass, stiffness and same-element
  face tensors;
* **numerical** -- a 15-point Gauss-Legendre rule along the face.  Used for
  the two tensors that couple *different* polynomial spaces
  (``INT_NODFUNC_TRI_FC_FC``, element x face) or *different elements*
  (``INT_NODFUNC_TRI_TRI_FC``, element x neighbour).

Array index order is kept identical to the Fortran so that the stage
comparison is a plain elementwise diff; only the 1-based/0-based convention
changes.  See ``docs/INDEXING.md``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .basis import Basis
from .constants import BC_PERIODIC
from .mesh import Mesh
from .quadrature import gauss_legendre_core, tri_quadrature

__all__ = ["Integrals", "build_integrals", "fac_div"]


# ---------------------------------------------------------------------------
def fac_div(ns: int, ne: int, ds: int, de: int) -> float:
    """``USD_Math.f90::FAC_DIV`` -- ``(ns..ne) / (ds..de)`` as a ratio of
    ascending products, evaluated term-by-term to avoid overflow."""
    numn = ne - ns + 1
    numd = de - ds + 1
    ans = 1.0
    for i in range(1, min(numn, numd) + 1):
        ans = ans * float(ns + i - 1) / float(ds + i - 1)
    if numn > numd:
        for i in range(numd + 1, numn + 1):
            ans = ans * float(ns + i - 1)
    if numd > numn:
        for i in range(numn + 1, numd + 1):
            ans = ans / float(ds + i - 1)
    return ans


def _affine_to_reference(x, y):
    """Coefficients of the physical -> reference map of a triangle.

    Written in the Fortran's exact form (``xi = B/A*xp + C/A*yp + D/A``).
    ``x``/``y`` are ``(..., 3)`` arrays of vertex coordinates.
    """
    x1, x2, x3 = x[..., 0], x[..., 1], x[..., 2]
    y1, y2, y3 = y[..., 0], y[..., 1], y[..., 2]
    A = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    B = y3 - y1
    C = x1 - x3
    D = (x3 - x1) * y1 - (y3 - y1) * x1
    E = (x3 - x1) * (y2 - y1) - (x2 - x1) * (y3 - y1)
    F = y2 - y1
    G = x1 - x2
    H = (x2 - x1) * y1 - (y2 - y1) * x1
    return (B / A, C / A, D / A), (F / E, G / E, H / E)


#: numpy's ``power`` special-cases small integral exponents and multiplies,
#: while gfortran's ``x**REAL(n,DBL)`` calls libm ``pow``.  The two differ by
#: one ulp on roughly 5% of inputs.  ``math.pow`` *is* libm ``pow``, so this
#: switch buys bit-identical face tensors at the cost of a Python-level loop
#: during setup only.  Turn it off for very large meshes.
EXACT_POW = True


def _powf(x, e):
    """``x ** float(e)`` with the Fortran's libm semantics."""
    import math

    if not EXACT_POW:
        return np.power(x, float(e))
    x = np.asarray(x, dtype=float)
    if e == 0:
        return np.ones_like(x)
    if e == 1:
        return x.copy()
    flat = x.ravel()
    out = np.fromiter((math.pow(v, float(e)) for v in flat),
                      dtype=np.float64, count=flat.size)
    return out.reshape(x.shape)


def _eval_mono(nodfun, mono, xi, eta):
    """``sum_k nodfun[m,k] * xi**a_k * eta**b_k`` for every ``m``.

    Accumulated in the Fortran's order (ascending monomial index, product
    associated as ``(c * xi**a) * eta**b``) rather than as a BLAS dot, so the
    result is bit-identical rather than merely equal to rounding.
    """
    shape = np.broadcast(xi, eta).shape
    out = np.zeros(shape + (nodfun.shape[0],))
    for k, (a, b) in enumerate(mono):
        xa = _powf(xi, a)
        yb = _powf(eta, b)
        out += (nodfun[:, k] * xa[..., None]) * yb[..., None]
    return out


def _contract2(C, E, order="il_outer"):
    """``ss[M, L] = sum C[L,IL] * C[M,IM] * E[IM,IL]`` in the Fortran's order.

    ``Integration.f90`` uses two different loop nestings for what is
    algebraically the same contraction; the summation order is the only
    difference and it shows up in the last bits, so both are offered.

    order='il_outer'  ->  DO IL; DO IM;  ss += C(L,IL)*C(M,IM)*E(IM,IL)
    order='im_outer'  ->  DO IM; DO IL;  ss += C(M,IM)*C(L,IL)*E(IM,IL)
    """
    n = C.shape[0]
    m = E.shape[0]
    out = np.zeros((n, n))
    if order == "il_outer":
        for il in range(m):
            for im in range(m):
                out += (C[:, il][None, :] * C[:, im][:, None]) * E[im, il]
    elif order == "im_outer":
        for im in range(m):
            for il in range(m):
                out += (C[:, im][:, None] * C[:, il][None, :]) * E[im, il]
    else:  # pragma: no cover
        raise ValueError(order)
    return out


def _contract1(C, e):
    """``ss[M] = sum_IM C[M,IM] * e[IM]``, accumulated in index order."""
    out = np.zeros(C.shape[0])
    for im in range(e.size):
        out += C[:, im] * e[im]
    return out


def _contract3(C, E3):
    """``ss[M,L,N] = sum C[N,IN]*C[L,IL]*C[M,IM]*E3[IM,IL,IN]``, IN/IL/IM."""
    n = C.shape[0]
    out = np.zeros((n, n, n))
    for in_ in range(n):
        for il in range(n):
            for im in range(n):
                out += ((C[:, in_][None, None, :] * C[:, il][None, :, None])
                        * C[:, im][:, None, None]) * E3[im, il, in_]
    return out


# ---------------------------------------------------------------------------
@dataclass
class Integrals:
    """All precomputed tensors, in Fortran index order.

    ``int_tri[L, I]``               ``INT_NODFUNC_TRI``
    ``int_tri_tri[M, L, I]``        ``INT_NODFUNC_TRI_TRI``
    ``int_tri_tri_x[I, M, L]``      ``INT_NODFUNC_TRI_TRI_X``
    ``int_tri_tri_y[I, M, L]``      ``INT_NODFUNC_TRI_TRI_Y``
    ``int_tri_fc[I, IL, M, L]``     ``INT_NODFUNC_TRI_FC``
    ``int_tri_fc_fc[M, L, I, IL]``  ``INT_NODFUNC_TRI_FC_FC``
    ``int_fc_fc[M, L, F]``          ``INT_NODFUNC_FC_FC``
    ``int_fc[M, F]``                ``INT_NODFUNC_FC``
    ``int_tri_tri_fc[I, IL, M, L]`` ``INT_NODFUNC_TRI_TRI_FC``
    """
    int_tri: np.ndarray
    int_tri_tri: np.ndarray
    int_tri_tri_x: np.ndarray
    int_tri_tri_y: np.ndarray
    int_tri_fc: np.ndarray
    int_tri_fc_fc: np.ndarray
    int_fc_fc: np.ndarray
    int_fc: np.ndarray
    int_tri_tri_fc: np.ndarray
    qua_abs_x: np.ndarray
    qua_abs_y: np.ndarray
    qua_wei: np.ndarray
    nodfun_qua_p: np.ndarray
    #: element-independent reference tensors, kept because the two big ones
    #: are just ``ref * 2 * area`` and materialising them per element wastes
    #: memory on refined grids
    ref_tri: np.ndarray
    ref_tri_tri: np.ndarray
    ref_tri_tri_tri: np.ndarray

    def int_tri_tri_tri(self, area: np.ndarray) -> np.ndarray:
        """``INT_NODFUNC_TRI_TRI_TRI(M,L,N,I)``; unused by the solver, kept
        for Gate-2 completeness."""
        return self.ref_tri_tri_tri[:, :, :, None] * (2.0 * area)[None, None, None, :]


def build_integrals(mesh: Mesh, basis: Basis, np_tri: int, np_fc: int,
                    bc_codes=None, bc_xoff=None, bc_yoff=None) -> Integrals:
    ndof_tri = basis.ndof_tri
    ndof_fc = basis.ndof_fc
    deg = basis.deg
    mono = basis.mono
    C = basis.nodfun_tri
    Cf = basis.nodfun_fc
    nodes = mesh.nodes
    n_tris = mesh.n_tris
    n_faces = mesh.n_faces
    area2 = 2.0 * mesh.tri_area

    ix = mono[:, 0]
    iy = mono[:, 1]

    # -- analytic monomial moment tables ----------------------------------
    e_vec = np.array([fac_div(1, iy[m], ix[m] + 1, ix[m] + iy[m] + 2)
                      for m in range(ndof_tri)])

    e_mass = np.empty((ndof_tri, ndof_tri))
    e_dx = np.zeros((ndof_tri, ndof_tri))
    e_dy = np.zeros((ndof_tri, ndof_tri))
    for m in range(ndof_tri):
        im, jm = int(ix[m]), int(iy[m])
        for l in range(ndof_tri):
            il, jl = int(ix[l]), int(iy[l])
            e_mass[m, l] = fac_div(1, jl + jm, il + im + 1, il + im + jl + jm + 2)
            if im > 0:
                e_dx[m, l] = im * fac_div(1, jl + jm, il + im, il + im - 1 + jl + jm + 2)
            if jm > 0:
                e_dy[m, l] = jm * fac_div(1, jl + jm - 1, il + im + 1,
                                          il + im - 1 + jl + jm + 2)

    ref_tri = _contract1(C, e_vec)                   # (ndof_tri,)
    ref_tri_tri = _contract2(C, e_mass, "il_outer")  # (M, L)
    ref_dx = _contract2(C, e_dx, "il_outer")
    ref_dy = _contract2(C, e_dy, "il_outer")

    # triple product (unused by the solver but part of the setup state)
    e3 = np.empty((ndof_tri, ndof_tri, ndof_tri))
    for m in range(ndof_tri):
        for l in range(ndof_tri):
            for n in range(ndof_tri):
                a = int(ix[m] + ix[l] + ix[n])
                b = int(iy[m] + iy[l] + iy[n])
                e3[m, l, n] = fac_div(1, b, a + 1, a + b + 2)
    ref_tri_tri_tri = _contract3(C, e3)

    int_tri = ref_tri[:, None] * area2[None, :]
    int_tri_tri = ref_tri_tri[:, :, None] * area2[None, None, :]

    # -- INT_NODFUNC_TRI_TRI_X / _Y ---------------------------------------
    x = nodes[mesh.tri_nodes, 0]
    y = nodes[mesh.tri_nodes, 1]
    x1, x2, x3 = x[:, 0], x[:, 1], x[:, 2]
    y1, y2, y3 = y[:, 0], y[:, 1], y[:, 2]
    ss = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)

    int_tri_tri_x = (((y3 - y1) / ss)[:, None, None] * ref_dx[None]
                     + ((y1 - y2) / ss)[:, None, None] * ref_dy[None]) * area2[:, None, None]
    int_tri_tri_y = (((x1 - x3) / ss)[:, None, None] * ref_dx[None]
                     + ((x2 - x1) / ss)[:, None, None] * ref_dy[None]) * area2[:, None, None]

    # -- INT_NODFUNC_TRI_FC (analytic edge integrals) ----------------------
    e_e1 = np.zeros((ndof_tri, ndof_tri))
    e_e2 = np.empty((ndof_tri, ndof_tri))
    e_e3 = np.zeros((ndof_tri, ndof_tri))
    sqrt2 = np.sqrt(2.0)
    for m in range(ndof_tri):
        im, jm = int(ix[m]), int(iy[m])
        for l in range(ndof_tri):
            il, jl = int(ix[l]), int(iy[l])
            if jm == 0 and jl == 0:
                e_e1[m, l] = 1.0 / float(im + il + 1)
            e_e2[m, l] = sqrt2 * fac_div(1, jm + jl, im + il + 1, im + il + jm + jl + 1)
            if im == 0 and il == 0:
                e_e3[m, l] = 1.0 / float(jm + jl + 1)

    r1 = _contract2(C, e_e1, "im_outer")
    r2 = _contract2(C, e_e2, "im_outer")
    r3 = _contract2(C, e_e3, "im_outer")

    len12 = np.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
    len23 = np.sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2))
    len31 = np.sqrt((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1))

    int_tri_fc = np.empty((n_tris, 3, ndof_tri, ndof_tri))
    int_tri_fc[:, 0] = len12[:, None, None] * r1[None]
    int_tri_fc[:, 1] = (len23 / sqrt2)[:, None, None] * r2[None]
    int_tri_fc[:, 2] = len31[:, None, None] * r3[None]

    # -- INT_NODFUNC_FC_FC and INT_NODFUNC_FC (analytic, on [-1,1]) --------
    e_line = np.zeros((ndof_fc, ndof_fc))
    for m in range(deg + 1):
        for l in range(deg + 1):
            if (m + l) % 2 == 0:
                e_line[m, l] = 2.0 / float(m + l + 1)
    ref_fc_fc = _contract2(Cf, e_line, "im_outer")
    ref_fc = _contract1(Cf, e_line[:, 0])
    half_len = mesh.face_len * 0.5
    int_fc_fc = ref_fc_fc[:, :, None] * half_len[None, None, :]
    int_fc = ref_fc[:, None] * half_len[None, :]

    # -- numerical face rules ----------------------------------------------
    yq, lpq = gauss_legendre_core(np_fc)
    n1, n2 = np_fc, np_fc + 1

    def face_rule(length):
        """GL rule on [0, length], written exactly as GaussLegendre does."""
        absc = (0.0 * (1.0 - yq) + length[:, None] * (1.0 + yq)) / 2.0
        weig = n2 * n2 * (length[:, None] - 0.0) / ((1.0 - yq * yq) * lpq * lpq) / n1 / n1
        return absc, weig

    int_tri_fc_fc = np.zeros((ndof_tri, ndof_fc, n_tris, 3))
    for j2 in range(3):
        fcid = mesh.tri_faces[:, j2]
        ln = mesh.face_len[fcid]
        absc, weig = face_rule(ln)
        fx1 = nodes[mesh.face_nodes[fcid, 0], 0][:, None]
        fy1 = nodes[mesh.face_nodes[fcid, 0], 1][:, None]
        fx2 = nodes[mesh.face_nodes[fcid, 1], 0][:, None]
        fy2 = nodes[mesh.face_nodes[fcid, 1], 1][:, None]
        lnc = ln[:, None]
        xp = (fx2 - fx1) / lnc * absc + fx1
        yp = (fy2 - fy1) / lnc * absc + fy1

        (bA, cA, dA), (fE, gE, hE) = _affine_to_reference(x, y)
        xi = bA[:, None] * xp + cA[:, None] * yp + dA[:, None]
        eta = fE[:, None] * xp + gE[:, None] * yp + hE[:, None]
        pm = _eval_mono(C, mono, xi, eta)             # (n_tris, np_fc, ndof_tri)

        t = 2.0 * absc / lnc - 1.0
        pl = np.zeros(t.shape + (ndof_fc,))           # (n_tris, np_fc, ndof_fc)
        for k in range(ndof_fc):
            pl += Cf[:, k] * _powf(t, k)[..., None]

        # accumulate over quadrature points in the Fortran's order
        acc = np.zeros((n_tris, ndof_tri, ndof_fc))
        for q in range(np_fc):
            acc += (pm[:, q, :, None] * pl[:, q, None, :]) * weig[:, q, None, None]
        int_tri_fc_fc[:, :, :, j2] = np.transpose(acc, (1, 2, 0))

    # -- INT_NODFUNC_TRI_TRI_FC --------------------------------------------
    bc_codes = np.asarray(bc_codes) if bc_codes is not None else np.zeros(0, dtype=int)
    bc_xoff = np.asarray(bc_xoff) if bc_xoff is not None else np.zeros(0)
    bc_yoff = np.asarray(bc_yoff) if bc_yoff is not None else np.zeros(0)

    int_tri_tri_fc = np.zeros((n_tris, 3, ndof_tri, ndof_tri))
    for i in range(n_tris):
        for j2 in range(3):
            fcid = int(mesh.tri_faces[i, j2])
            bc = int(mesh.face_bc[fcid])
            ln = float(mesh.face_len[fcid])
            fx1 = nodes[mesh.face_nodes[fcid, 0], 0]
            fy1 = nodes[mesh.face_nodes[fcid, 0], 1]
            fx2 = nodes[mesh.face_nodes[fcid, 1], 0]
            fy2 = nodes[mesh.face_nodes[fcid, 1], 1]
            absc, weig = face_rule(np.array([ln]))
            absc = absc[0]
            weig = weig[0]
            xp = (fx2 - fx1) / ln * absc + fx1
            yp = (fy2 - fy1) / ln * absc + fy1

            xs = nodes[mesh.tri_nodes[i], 0]
            ys = nodes[mesh.tri_nodes[i], 1]
            (bA, cA, dA), (fE, gE, hE) = _affine_to_reference(xs, ys)
            pm = _eval_mono(C, mono, bA * xp + cA * yp + dA, fE * xp + gE * yp + hE)

            periodic = bc >= 0 and bc_codes.size and bc_codes[bc] == BC_PERIODIC
            if bc < 0 or periodic:
                if bc < 0:
                    t1, t2 = mesh.face_tri[fcid]
                    other = t2 if t1 == i else t1
                    ox = oy = 0.0
                else:
                    pair = int(mesh.face_pair[fcid])
                    bc_p = int(mesh.face_bc[pair])
                    other, _ = mesh.adjacent(pair)
                    ox, oy = float(bc_xoff[bc_p]), float(bc_yoff[bc_p])
                xo = nodes[mesh.tri_nodes[other], 0] + ox
                yo = nodes[mesh.tri_nodes[other], 1] + oy
                (bA2, cA2, dA2), (fE2, gE2, hE2) = _affine_to_reference(xo, yo)
                pl = _eval_mono(C, mono, bA2 * xp + cA2 * yp + dA2,
                                fE2 * xp + gE2 * yp + hE2)
                acc = np.zeros((ndof_tri, ndof_tri))
                for q in range(np_fc):
                    acc += (pm[q][:, None] * pl[q][None, :]) * weig[q]
                int_tri_tri_fc[i, j2] = acc
            else:
                acc = np.zeros(ndof_tri)
                for q in range(np_fc):
                    acc += pm[q] * weig[q]
                int_tri_tri_fc[i, j2] = acc[:, None]

    # -- element quadrature ------------------------------------------------
    qax, qay, qw = tri_quadrature(np_tri)
    nodfun_qua_p = _eval_mono(C, mono, qax, qay).T    # (ndof_tri, np_tri)

    return Integrals(
        int_tri=int_tri, int_tri_tri=int_tri_tri,
        int_tri_tri_x=int_tri_tri_x, int_tri_tri_y=int_tri_tri_y,
        int_tri_fc=int_tri_fc, int_tri_fc_fc=int_tri_fc_fc,
        int_fc_fc=int_fc_fc, int_fc=int_fc, int_tri_tri_fc=int_tri_tri_fc,
        qua_abs_x=qax, qua_abs_y=qay, qua_wei=qw, nodfun_qua_p=nodfun_qua_p,
        ref_tri=ref_tri, ref_tri_tri=ref_tri_tri, ref_tri_tri_tri=ref_tri_tri_tri,
    )
