"""§6.1 unit level: the precomputed integral tensors, checked against
independently computed quadrature rather than against the Fortran."""
from __future__ import annotations

import numpy as np
import pytest

from pybte.basis import eval_tri_basis
from pybte.integration import build_integrals, fac_div
from pybte.quadrature import tri_quadrature


@pytest.fixture(scope="module")
def integ(shipped_mesh, basis3):
    return build_integrals(shipped_mesh, basis3, 12, 15,
                           bc_codes=[1, 1, 1, 1], bc_xoff=[0] * 4, bc_yoff=[0] * 4)


def test_fac_div():
    from math import factorial
    for ns, ne, ds, de in [(1, 0, 1, 2), (1, 3, 2, 6), (1, 5, 3, 4), (2, 2, 1, 9)]:
        num = np.prod([float(k) for k in range(ns, ne + 1)]) if ne >= ns else 1.0
        den = np.prod([float(k) for k in range(ds, de + 1)]) if de >= ds else 1.0
        assert fac_div(ns, ne, ds, de) == pytest.approx(num / den, rel=1e-13)
    # the shape actually used: a! b! / (a+b+2)!
    for a in range(5):
        for b in range(5):
            want = factorial(a) * factorial(b) / factorial(a + b + 2)
            assert fac_div(1, b, a + 1, a + b + 2) == pytest.approx(want, rel=1e-13)


def test_mass_matrix_row_sums_equal_area(integ, shipped_mesh):
    """sum_L int phi_M phi_L = int phi_M, and summing over M gives the area."""
    rowsum = integ.int_tri_tri.sum(axis=1)          # (M, I)
    assert np.allclose(rowsum, integ.int_tri, rtol=1e-12)
    assert np.allclose(integ.int_tri.sum(axis=0), shipped_mesh.tri_area, rtol=1e-12)


def test_mass_matrix_symmetric(integ):
    """Symmetric to 1e-13 relative, not exactly: the Fortran contracts the
    monomial moments with the nodal coefficients in a fixed loop order, so
    entry (M,L) and entry (L,M) are the same terms summed differently.  The
    sweep indexes the transpose (``INT_NODFUNC_TRI_TRI(L,M,TRID)``) in the
    source and the un-transposed form in the operator, so the port has to
    keep both index orders rather than assume symmetry."""
    m = integ.int_tri_tri
    scale = np.abs(m).max()
    assert np.abs(m - np.transpose(m, (1, 0, 2))).max() <= 1e-13 * scale


def test_mass_matrix_against_quadrature(integ, shipped_mesh, basis3):
    """Independent check with a 16-point triangle rule."""
    x, y, w = tri_quadrature(16)
    phi = eval_tri_basis(basis3.nodfun_tri, basis3.mono, x, y)   # (np, ndof)
    # the reference triangle has area 1/2 and the rule's weights sum to 1
    ref = 0.5 * np.einsum("qm,ql,q->ml", phi, phi, w, optimize=True)
    for i in [0, 37, 199]:
        got = integ.int_tri_tri[:, :, i]
        want = ref * 2.0 * shipped_mesh.tri_area[i]
        assert np.abs(got - want).max() <= 1e-12 * np.abs(want).max()


def test_stiffness_partition(integ):
    """sum_L int (d phi_M/dx) phi_L = int d phi_M/dx, and summing over M gives
    zero because sum_M phi_M = 1 has zero derivative."""
    assert np.allclose(integ.int_tri_tri_x.sum(axis=1), 0.0, atol=1e-12)
    assert np.allclose(integ.int_tri_tri_y.sum(axis=1), 0.0, atol=1e-12)


def test_divergence_theorem(integ, shipped_mesh):
    """int_T d phi_M/dx dT = oint phi_M n_x ds, elementwise.

    This ties the volume and the face tensors together and would catch a
    wrong edge length, a flipped normal or a mis-scaled reference map.
    """
    vol_x = integ.int_tri_tri_x.sum(axis=2)           # (I, M): int dphi_M/dx
    vol_y = integ.int_tri_tri_y.sum(axis=2)
    surf_x = np.zeros_like(vol_x)
    surf_y = np.zeros_like(vol_y)
    for il in range(3):
        edge = integ.int_tri_fc[:, il].sum(axis=2)    # (I, M): oint phi_M ds
        surf_x += edge * shipped_mesh.tri_normal[:, il, 0][:, None]
        surf_y += edge * shipped_mesh.tri_normal[:, il, 1][:, None]
    assert np.allclose(vol_x, surf_x, atol=1e-13)
    assert np.allclose(vol_y, surf_y, atol=1e-13)


def test_face_tensors_partition(integ, shipped_mesh):
    """oint phi_M phi_L summed over both indices is the edge length."""
    for il in range(3):
        tot = integ.int_tri_fc[:, il].sum(axis=(1, 2))
        want = shipped_mesh.face_len[shipped_mesh.tri_faces[:, il]]
        assert np.allclose(tot, want, rtol=1e-12)


def test_tri_fc_fc_partition(integ, shipped_mesh):
    """The element x face tensor also sums to the edge length."""
    for il in range(3):
        tot = integ.int_tri_fc_fc[:, :, :, il].sum(axis=(0, 1))
        want = shipped_mesh.face_len[shipped_mesh.tri_faces[:, il]]
        assert np.allclose(tot, want, rtol=1e-12)


def test_neighbour_face_tensor_is_consistent(integ, shipped_mesh):
    """On an interior face the element x neighbour tensor must reduce to the
    same edge length, because both bases are partitions of unity there."""
    m = shipped_mesh
    for i in range(m.n_tris):
        for il in range(3):
            if m.face_bc[m.tri_faces[i, il]] >= 0:
                continue
            tot = integ.int_tri_tri_fc[i, il].sum()
            assert tot == pytest.approx(m.face_len[m.tri_faces[i, il]], rel=1e-12)


def test_fc_fc_and_fc(integ, shipped_mesh):
    assert np.allclose(integ.int_fc_fc.sum(axis=(0, 1)), shipped_mesh.face_len, rtol=1e-12)
    assert np.allclose(integ.int_fc.sum(axis=0), shipped_mesh.face_len, rtol=1e-12)
    f = integ.int_fc_fc
    assert np.allclose(f, np.transpose(f, (1, 0, 2)), rtol=1e-13, atol=1e-18)
