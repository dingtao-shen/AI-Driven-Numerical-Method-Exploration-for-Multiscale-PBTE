"""Nodal (Lagrange) bases -- the port of ``Basis_Function.f90``.

Both bases are stored as **monomial coefficients**, exactly as the Fortran
does, because every integral downstream is assembled by contracting those
coefficients with an analytic monomial-moment table.

Face basis, on the reference segment ``[-1, 1]``::

    node_fc[i]        = -1 + 2*i/(NDOF_FC-1)
    phi_i(t)          = sum_k nodfun_fc[i, k] * t**k

Element basis, on the reference triangle (0,0)-(1,0)-(0,1)::

    node_tri[c]       = (i/DEG, j/DEG) for j = 0..DEG, i = 0..DEG-j
    phi_c(xi, eta)    = sum_k nodfun_tri[c, k] * xi**ix[k] * eta**iy[k]

with the monomial ordering fixed by :func:`monomial_exponents`: grouped by
total degree, and within a degree by the eta exponent.  Index of
``xi**a eta**b`` is ``(a+b)*(a+b+1)//2 + b``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["Basis", "build_basis", "monomial_exponents", "monomial_index",
           "basis_function_fc", "basis_function_tri", "eval_tri_basis"]


def monomial_exponents(deg: int) -> np.ndarray:
    """``(ndof, 2)`` array of ``(xi_exp, eta_exp)`` in the Fortran's order."""
    out = []
    for k in range(deg + 1):
        for j in range(k + 1):
            out.append((k - j, j))
    return np.asarray(out, dtype=np.int64)


def monomial_index(a: int, b: int) -> int:
    """0-based index of ``xi**a eta**b``."""
    return (a + b) * (a + b + 1) // 2 + b


# ---------------------------------------------------------------------------
def basis_function_fc(nn: int):
    """``(node_ref, nodfun_ref)`` for the 1-D nodal basis on ``[-1, 1]``.

    Transcribes ``Basis_Function_FC``: the Lagrange product is expanded by
    repeatedly multiplying by ``b0 + b1*t`` with
    ``b0 = (2j-nn-1)/(2(j-i))`` and ``b1 = -(nn-1)/(2(j-i))``.
    """
    if nn < 2:
        raise ValueError("NDOF_FC = DEG+1 must be >= 2")
    node = np.array([-1.0 + 2.0 * i / (nn - 1) for i in range(nn)])
    coef = np.zeros((nn, nn))
    for i in range(nn):                      # i is 0-based, Fortran's I-1
        coef[i, 0] = 1.0
        k = 0
        for j in range(nn):
            if j == i:
                continue
            a = coef[i, :k + 1].copy()
            b0 = float(2 * (j + 1) - nn - 1) / float(2 * (j - i))
            b1 = -float(nn - 1) / float(2 * (j - i))
            coef[i, 0] = a[0] * b0
            coef[i, k + 1] = a[k] * b1
            for ell in range(1, k + 1):
                coef[i, ell] = a[ell] * b0 + a[ell - 1] * b1
            k += 1
    return node, coef


def basis_function_tri(rr: int):
    """``(node_ref, nodfun_ref)`` for the nodal basis on the reference triangle.

    Transcribes ``Basis_Function_TRI``.  The shape function of node
    ``(i/rr, j/rr)`` is the triple product

        prod_{l<i} (rr*xi - l)/(i - l)
      * prod_{m<j} (rr*eta - m)/(j - m)
      * prod_{n<k} (rr*(1-xi-eta) - n)/(k - n)      with k = rr - i - j
    """
    nn = (rr + 1) * (rr + 2) // 2
    node = np.zeros((nn, 2))
    coef = np.zeros((nn, nn))

    c = 0
    for j in range(rr + 1):
        for i in range(rr - j + 1):
            node[c, 0] = float(i) / float(rr)
            node[c, 1] = float(j) / float(rr)
            c += 1

    c = 0
    for j in range(rr + 1):
        for i in range(rr - j + 1):
            k = rr - i - j
            a = np.zeros(nn)
            b = np.zeros(nn)

            # --- prod_{l<i} (rr*xi - l)/(i - l), stored at index of xi**l ---
            if i != 0:
                a[1] = float(rr) / float(i)          # a(2) in Fortran: xi**1
                for ell in range(1, i):              # L = 1 .. i-1
                    tmp = a.copy()
                    d0 = float(ell) / float(ell - i)
                    d1 = -float(rr) / float(ell - i)
                    a[0] = tmp[0] * d0
                    a[(ell + 1) * (ell + 2) // 2] = tmp[ell * (ell + 1) // 2] * d1
                    for m in range(1, ell + 1):
                        a[m * (m + 1) // 2] = (tmp[m * (m + 1) // 2] * d0
                                               + tmp[(m - 1) * m // 2] * d1)
            else:
                a[0] = 1.0

            # --- prod_{m<j} (rr*eta - m)/(j - m), stored at index of eta**m --
            if j != 0:
                b[2] = float(rr) / float(j)          # b(3) in Fortran: eta**1
                for m in range(1, j):                # M = 1 .. j-1
                    tmp = b.copy()
                    d0 = float(m) / float(m - j)
                    d1 = -float(rr) / float(m - j)
                    b[0] = tmp[0] * d0
                    b[(m + 2) * (m + 3) // 2 - 1] = tmp[(m + 1) * (m + 2) // 2 - 1] * d1
                    for ell in range(1, m + 1):
                        b[(ell + 1) * (ell + 2) // 2 - 1] = (
                            tmp[(ell + 1) * (ell + 2) // 2 - 1] * d0
                            + tmp[ell * (ell + 1) // 2 - 1] * d1)
            else:
                b[0] = 1.0

            # --- the product a(xi) * b(eta) --------------------------------
            for ell in range(i + 1):
                for m in range(j + 1):
                    pos = (ell + m) * (ell + m + 1) // 2 + m
                    coef[c, pos] += (a[ell * (ell + 1) // 2]
                                     * b[(m + 1) * (m + 2) // 2 - 1])

            # --- multiply by prod_{n<k} (rr*(1-xi-eta) - n)/(k - n) ---------
            if k != 0:
                for n in range(k):
                    tmp = coef[c, :].copy()
                    d0 = float(n - rr) / float(n - k)
                    d1 = float(rr) / float(n - k)
                    d2 = float(rr) / float(n - k)
                    for cc in range(rr + 1):
                        for m in range(cc + 1):
                            ell = cc - m
                            pos = (ell + m) * (ell + m + 1) // 2 + m
                            coef[c, pos] = tmp[pos] * d0
                            if ell > 0:
                                pos1 = (ell - 1 + m) * (ell + m) // 2 + m
                                coef[c, pos] += tmp[pos1] * d1
                            if m > 0:
                                pos2 = (ell + m - 1) * (ell + m) // 2 + m - 1
                                coef[c, pos] += tmp[pos2] * d2
            c += 1

    return node, coef


# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Basis:
    deg: int
    node_fc: np.ndarray        # (ndof_fc,)
    nodfun_fc: np.ndarray      # (ndof_fc, ndof_fc)
    node_tri: np.ndarray       # (ndof_tri, 2)
    nodfun_tri: np.ndarray     # (ndof_tri, ndof_tri)
    mono: np.ndarray           # (ndof_tri, 2) monomial exponents

    @property
    def ndof_tri(self) -> int:
        return self.nodfun_tri.shape[0]

    @property
    def ndof_fc(self) -> int:
        return self.nodfun_fc.shape[0]

    def eval_tri(self, xi, eta) -> np.ndarray:
        return eval_tri_basis(self.nodfun_tri, self.mono, xi, eta)

    def eval_fc(self, t) -> np.ndarray:
        t = np.atleast_1d(np.asarray(t, dtype=float))
        powers = np.stack([np.power(t, float(k)) for k in range(self.ndof_fc)], axis=-1)
        return powers @ self.nodfun_fc.T


def build_basis(deg: int) -> Basis:
    node_fc, nodfun_fc = basis_function_fc(deg + 1)
    node_tri, nodfun_tri = basis_function_tri(deg)
    return Basis(deg=deg, node_fc=node_fc, nodfun_fc=nodfun_fc,
                 node_tri=node_tri, nodfun_tri=nodfun_tri,
                 mono=monomial_exponents(deg))


def eval_tri_basis(nodfun_tri: np.ndarray, mono: np.ndarray, xi, eta) -> np.ndarray:
    """Evaluate every element shape function at one or many points.

    The Fortran writes ``xi**REAL(IL,DBL)``, i.e. C ``pow`` with a
    floating-point exponent, so we do the same rather than repeated
    multiplication -- the two differ in the last bit for degree >= 3.
    """
    xi = np.asarray(xi, dtype=float)
    eta = np.asarray(eta, dtype=float)
    scalar = xi.ndim == 0
    xi = np.atleast_1d(xi)
    eta = np.atleast_1d(eta)
    terms = np.empty((xi.size, mono.shape[0]))
    for k, (a, b) in enumerate(mono):
        terms[:, k] = np.power(xi, float(a)) * np.power(eta, float(b))
    out = terms @ nodfun_tri.T
    return out[0] if scalar else out
