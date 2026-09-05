"""§6.1 unit level: mesh topology and geometry."""
from __future__ import annotations

import numpy as np
import pytest

from pybte.mesh import build_mesh
from pybte.mesh.gmsh_reader import read_gmsh22

from ..conftest import MESH


def test_counts_and_identity(shipped_mesh):
    m = shipped_mesh
    assert m.n_nodes == 121
    assert m.n_tris == 200
    assert m.n_faces == 320
    assert m.n_faces_b == 40
    assert m.n_faces_i == 280
    # 3*T = 2*I + B for a conforming triangulation
    assert 3 * m.n_tris == 2 * m.n_faces_i + m.n_faces_b


def test_every_interior_face_has_two_neighbours(shipped_mesh):
    m = shipped_mesh
    interior = m.face_bc < 0
    assert np.all(m.face_tri[interior] >= 0)
    assert np.all(m.face_lfc[interior] >= 0)
    # ... and exactly one for a boundary face
    bnd = ~interior
    assert np.all((m.face_tri[bnd] >= 0).sum(axis=1) == 1)


def test_boundary_faces_carry_a_valid_bcid(shipped_mesh):
    m = shipped_mesh
    bnd = m.face_bc >= 0
    assert np.all(m.face_bc[bnd] < 4)
    assert int(bnd.sum()) == 40
    # 10 faces per wall on an 11x11 node grid
    for k in range(4):
        assert int(np.count_nonzero(m.face_bc == k)) == 10


def test_boundary_faces_come_first(shipped_mesh):
    """``Calculate_FLUX_WALL`` loops ``1..N_FCS_B`` and indexes FLUX_WALL by
    face id, which is only valid if the boundary faces are the first block.
    Our own code does not rely on it, but a mesh that breaks the assumption
    would silently mean something different to the reference."""
    m = shipped_mesh
    assert np.all(m.face_bc[:m.n_faces_b] >= 0)
    assert np.all(m.face_bc[m.n_faces_b:] < 0)


def test_areas_and_normals(shipped_mesh):
    m = shipped_mesh
    assert m.tri_area.sum() == pytest.approx(1.0, rel=1e-13)
    assert np.all(m.tri_area > 0)
    n = m.tri_normal
    assert np.allclose(np.sum(n ** 2, axis=2), 1.0, atol=1e-14)

    # normals point away from the centroid
    cen = m.nodes[m.tri_nodes].mean(axis=1)
    for j in range(3):
        a = m.nodes[m.tri_nodes[:, j]]
        b = m.nodes[m.tri_nodes[:, (j + 1) % 3]]
        mid = 0.5 * (a + b)
        assert np.all(np.sum(n[:, j] * (mid - cen), axis=1) > 0)


def test_face_length_matches_edges(shipped_mesh):
    m = shipped_mesh
    for i in range(m.n_tris):
        for j in range(3):
            f = m.tri_faces[i, j]
            a = m.nodes[m.tri_nodes[i, j]]
            b = m.nodes[m.tri_nodes[i, (j + 1) % 3]]
            assert m.face_len[f] == pytest.approx(np.linalg.norm(a - b), rel=1e-14)


def test_neighbour_symmetry(shipped_mesh):
    m = shipped_mesh
    for i in range(m.n_tris):
        for j in range(3):
            nb = m.tri_neighbour[i, j]
            if nb >= 0:
                assert i in m.tri_neighbour[nb]


def test_hmin(shipped_mesh):
    m = shipped_mesh
    assert m.hmin == pytest.approx(m.tri_hmin.min(), rel=1e-15)
    assert 0 < m.hmin < 1


def test_reader_rejects_unsupported_element_type(tmp_path):
    """The Fortran BACKSPACEs and never re-reads for types other than 1 and 2,
    desynchronising the whole element list.  We refuse instead."""
    text = MESH.read_text().replace("    41 2 2 15 10", "    41 15 2 15 10", 1)
    p = tmp_path / "bad.msh"
    p.write_text(text)
    with pytest.raises(ValueError, match="unsupported type"):
        read_gmsh22(p)


def test_untagged_boundary_is_rejected(tmp_path):
    """A boundary face whose physical id is not in the BC list would end up
    with BCID 0 (= interior) and produce a one-sided interior face."""
    with pytest.raises(ValueError):
        build_mesh(MESH, [11, 12, 13], [1, 1, 1], ["S", "N", "E"])
