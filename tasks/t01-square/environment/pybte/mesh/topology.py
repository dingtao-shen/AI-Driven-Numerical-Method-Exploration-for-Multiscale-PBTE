"""Face extraction, connectivity and periodic pairing.

This is the bookkeeping half of ``Spatial_Mesh.f90::Init_Spatial_Grid``.  The
*ordering* produced here is load-bearing: triangle ids feed ``TRI_ORDER``'s
tie-breaking, so it must match the Fortran element-for-element.

Face numbering, exactly as the Fortran produces it:

1. elements are visited in gmsh file order;
2. a type-1 (line) element always appends a new face, tagged with the index of
   the first matching ``BC_PHYID``;
3. a type-2 (triangle) element appends a face for each of its three edges that
   is not already present as an undirected node pair.

Indexing convention (docs/INDEXING.md): everything here is **0-based**, and
``-1`` means "absent" where the Fortran uses 0.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import BC_SYMMETRY, BC_THERMALISING
from .geometry import triangle_geometry
from .gmsh_reader import LINE_ELEMENT, TRIANGLE_ELEMENT, GmshMesh, read_gmsh22

__all__ = ["Mesh", "build_mesh"]


@dataclass
class Mesh:
    """The spatial discretisation.

    Attributes
    ----------
    nodes : (n_nodes, 2)
    tri_nodes : (n_tris, 3) int   0-based node ids, gmsh order
    tri_faces : (n_tris, 3) int   0-based face id of local edges 1,2,3
    tri_area : (n_tris,)
    tri_normal : (n_tris, 3, 2)   outward unit normal per local edge
    tri_hmin : (n_tris,)          smallest height 2A/L
    hmin : float                  ``min(1.0, min over all heights)``
    face_nodes : (n_faces, 2) int
    face_bc : (n_faces,) int      0-based BC index, -1 for interior
    face_pair : (n_faces,) int    periodic partner, or the wall/symmetry queue
                                  index for BC types 1 and 4, else -1
    face_tri : (n_faces, 2) int   [triangle(+), triangle(-)], -1 if absent
    face_lfc : (n_faces, 2) int   local edge id (0-based) in triangle(+)/(-)
    face_len : (n_faces,)
    tri_neighbour : (n_tris, 3) int  neighbour across each local edge, -1 if
                                  the edge is on a boundary
    """
    nodes: np.ndarray
    tri_nodes: np.ndarray
    tri_faces: np.ndarray
    tri_area: np.ndarray
    tri_normal: np.ndarray
    tri_hmin: np.ndarray
    hmin: float
    face_nodes: np.ndarray
    face_bc: np.ndarray
    face_pair: np.ndarray
    face_tri: np.ndarray
    face_lfc: np.ndarray
    face_len: np.ndarray
    tri_neighbour: np.ndarray
    bc_codes: np.ndarray

    # -- sizes ------------------------------------------------------------
    @property
    def n_nodes(self) -> int:
        return self.nodes.shape[0]

    @property
    def n_tris(self) -> int:
        return self.tri_nodes.shape[0]

    @property
    def n_faces(self) -> int:
        return self.face_nodes.shape[0]

    @property
    def n_faces_b(self) -> int:
        return int(np.count_nonzero(self.face_bc >= 0))

    @property
    def n_faces_i(self) -> int:
        return self.n_faces - self.n_faces_b

    # -- Fortran-shaped views, for cross-validation only -------------------
    def fortran_triangles_tag(self) -> np.ndarray:
        out = np.zeros((self.n_tris, 6), dtype=np.int32)
        out[:, 0:3] = self.tri_nodes + 1
        out[:, 3:6] = self.tri_faces + 1
        return out

    def fortran_triangles_inf(self) -> np.ndarray:
        out = np.zeros((self.n_tris, 7), dtype=np.float64)
        out[:, 0] = self.tri_area
        out[:, 1:7] = self.tri_normal.reshape(self.n_tris, 6)
        return out

    def fortran_faces_tag(self) -> np.ndarray:
        out = np.zeros((self.n_faces, 8), dtype=np.int32)
        out[:, 0:2] = self.face_nodes + 1
        out[:, 2] = self.face_bc + 1
        out[:, 3] = self.face_pair + 1
        out[:, 4] = self.face_tri[:, 0] + 1
        out[:, 5] = self.face_tri[:, 1] + 1
        out[:, 6] = self.face_lfc[:, 0] + 1
        out[:, 7] = self.face_lfc[:, 1] + 1
        return out

    # -- helpers -----------------------------------------------------------
    def adjacent(self, face: int):
        """``(triangle, local_edge)`` of the single element touching a
        boundary face, using the Fortran's plus-then-minus preference."""
        if self.face_tri[face, 0] < 0:
            return int(self.face_tri[face, 1]), int(self.face_lfc[face, 1])
        return int(self.face_tri[face, 0]), int(self.face_lfc[face, 0])


def build_mesh(path, bc_phyid, bc_codes, bc_names=None,
               bc_xoff=None, bc_yoff=None) -> Mesh:
    """Read a gmsh 2.2 file and build the full topology.

    Parameters
    ----------
    path : str or Path
    bc_phyid : sequence[int]     physical ids, in ``control.in`` order
    bc_codes : sequence[int]     ``BC_TYP`` codes, same order
    bc_names : sequence[str]     used only to find the periodic Master/Slave
    bc_xoff, bc_yoff : sequence[float]
    """
    g: GmshMesh = read_gmsh22(path)
    nbc = len(bc_phyid)
    bc_phyid = np.asarray(bc_phyid, dtype=np.int64)
    bc_codes = np.asarray(bc_codes, dtype=np.int64)
    bc_names = list(bc_names) if bc_names is not None else [""] * nbc
    bc_xoff = np.asarray(bc_xoff if bc_xoff is not None else np.zeros(nbc), dtype=float)
    bc_yoff = np.asarray(bc_yoff if bc_yoff is not None else np.zeros(nbc), dtype=float)

    face_nodes: list[tuple[int, int]] = []
    face_bc: list[int] = []
    key_to_face: dict[tuple[int, int], int] = {}
    tri_nodes: list[tuple[int, int, int]] = []

    for etype, phyid, nds in zip(g.elem_type, g.elem_phyid, g.elem_nodes):
        if etype == LINE_ELEMENT:
            n1, n2 = nds
            bc = -1
            for j in range(nbc):
                if phyid == bc_phyid[j]:
                    bc = j          # Fortran uses CYCLE, so the last match wins
            key = (min(n1, n2), max(n1, n2))
            if key in key_to_face:
                raise ValueError(
                    f"{path}: duplicate boundary edge {(n1 + 1, n2 + 1)}; the "
                    "reference reader would create two faces for it")
            key_to_face[key] = len(face_nodes)
            face_nodes.append((n1, n2))
            face_bc.append(bc)
        elif etype == TRIANGLE_ELEMENT:
            nd = list(nds) + [nds[0]]
            tri_nodes.append(tuple(nds))
            for j in range(3):
                a, b = nd[j], nd[j + 1]
                key = (min(a, b), max(a, b))
                if key not in key_to_face:
                    key_to_face[key] = len(face_nodes)
                    face_nodes.append((a, b))
                    face_bc.append(-1)

    n_tris = len(tri_nodes)
    n_faces = len(face_nodes)
    if n_tris == 0:
        raise ValueError(f"{path}: no triangles")

    tri_nodes_a = np.asarray(tri_nodes, dtype=np.int32)
    face_nodes_a = np.asarray(face_nodes, dtype=np.int32)
    face_bc_a = np.asarray(face_bc, dtype=np.int32)
    face_pair = np.full(n_faces, -1, dtype=np.int32)

    # ---- periodic pairing (BC_NAME 'Master' / 'Slave') -------------------
    id_master = next((i for i, s in enumerate(bc_names) if str(s).strip() == "Master"), -1)
    id_slave = next((i for i, s in enumerate(bc_names) if str(s).strip() == "Slave"), -1)
    if id_master >= 0:
        if id_slave < 0:
            raise ValueError("a 'Master' boundary needs a matching 'Slave'")
        masters = np.flatnonzero(face_bc_a == id_master)
        slaves = np.flatnonzero(face_bc_a == id_slave)
        xoff, yoff = bc_xoff[id_master], bc_yoff[id_master]
        for m in masters:
            a, b = face_nodes_a[m]
            pa = (g.nodes[a, 0] + xoff, g.nodes[a, 1] + yoff)
            pb = (g.nodes[b, 0] + xoff, g.nodes[b, 1] + yoff)
            for s in slaves:
                c, d = face_nodes_a[s]
                pc = (g.nodes[c, 0], g.nodes[c, 1])
                pd = (g.nodes[d, 0], g.nodes[d, 1])
                if (pa == pc and pb == pd) or (pa == pd and pb == pc):
                    face_pair[m] = s
                    face_pair[s] = m
                    break
        unpaired = [int(i) for i in np.concatenate([masters, slaves]) if face_pair[i] < 0]
        if unpaired:
            raise ValueError(f"unpaired periodic faces: {unpaired[:8]}")

    # ---- geometry --------------------------------------------------------
    area, normal, hmin_tri, hmin = triangle_geometry(g.nodes, tri_nodes_a)

    fn = g.nodes[face_nodes_a]                              # (n_faces, 2, 2)
    d = fn[:, 0, :] - fn[:, 1, :]
    face_len = np.sqrt(d[:, 0] * d[:, 0] + d[:, 1] * d[:, 1])

    # ---- connect faces and triangle edges --------------------------------
    tri_faces = np.full((n_tris, 3), -1, dtype=np.int32)
    face_tri = np.full((n_faces, 2), -1, dtype=np.int32)
    face_lfc = np.full((n_faces, 2), -1, dtype=np.int32)

    for i in range(n_tris):
        nd = list(tri_nodes_a[i]) + [tri_nodes_a[i][0]]
        for j in range(3):
            a, b = nd[j], nd[j + 1]
            k = key_to_face[(min(a, b), max(a, b))]
            tri_faces[i, j] = k
            if face_nodes_a[k, 0] == a and face_nodes_a[k, 1] == b:
                face_tri[k, 0] = i          # triangle(+)
                face_lfc[k, 0] = j
            else:
                face_tri[k, 1] = i          # triangle(-)
                face_lfc[k, 1] = j

    # ---- wall / symmetry queue (overwrites FACES_TAG(:,4) as in Fortran) --
    counter = 0
    for k in range(n_faces):
        bc = face_bc_a[k]
        if bc >= 0 and bc_codes[bc] in (BC_THERMALISING, BC_SYMMETRY):
            face_pair[k] = counter
            counter += 1

    # ---- neighbours across each local edge --------------------------------
    tri_neighbour = np.full((n_tris, 3), -1, dtype=np.int32)
    for i in range(n_tris):
        for j in range(3):
            k = tri_faces[i, j]
            t1, t2 = face_tri[k]
            if i == t1:
                tri_neighbour[i, j] = t2
            elif i == t2:
                tri_neighbour[i, j] = t1

    mesh = Mesh(nodes=g.nodes, tri_nodes=tri_nodes_a, tri_faces=tri_faces,
                tri_area=area, tri_normal=normal, tri_hmin=hmin_tri, hmin=hmin,
                face_nodes=face_nodes_a, face_bc=face_bc_a, face_pair=face_pair,
                face_tri=face_tri, face_lfc=face_lfc, face_len=face_len,
                tri_neighbour=tri_neighbour, bc_codes=bc_codes)

    # sanity: an interior face must have both neighbours, a boundary face one
    interior = mesh.face_bc < 0
    if np.any((mesh.face_tri[interior] < 0).any(axis=1)):
        bad = np.flatnonzero(interior & (mesh.face_tri < 0).any(axis=1))
        raise ValueError(f"interior faces with a missing neighbour: {bad[:8].tolist()}. "
                         "Either the mesh has a hole or a boundary is untagged.")
    boundary = ~interior
    if np.any((mesh.face_tri[boundary] >= 0).all(axis=1)):
        bad = np.flatnonzero(boundary & (mesh.face_tri >= 0).all(axis=1))
        raise ValueError(f"boundary faces with two neighbours: {bad[:8].tolist()}")
    if np.any(mesh.face_bc[boundary] < 0):
        raise ValueError("a boundary face carries no BC index")

    return mesh
