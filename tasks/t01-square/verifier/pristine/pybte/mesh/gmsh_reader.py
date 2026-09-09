"""Gmsh 2.2 ASCII reader -- the parsing half of ``Spatial_Mesh.f90``.

The Fortran reads only ``$Nodes`` and ``$Elements``, keeps element types 1
(2-node line) and 2 (3-node triangle), and assumes exactly two integer tags
per element (physical id first, geometrical id second).  Anything else in the
file is ignored.  We reproduce that, but *fail loudly* on the cases where the
Fortran would silently misbehave:

* an element of a type other than 1 or 2 leaves the Fortran's file pointer on
  the wrong line (it BACKSPACEs and never re-reads), corrupting everything
  after it -- we raise instead;
* fewer than two tags would make the Fortran read node ids as tags.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["GmshMesh", "read_gmsh22"]

LINE_ELEMENT = 1
TRIANGLE_ELEMENT = 2


@dataclass
class GmshMesh:
    """Raw contents of the file, in file order."""
    nodes: np.ndarray            # (n_nodes, 2)  x, y  (z discarded)
    elem_type: np.ndarray        # (n_elements,) int
    elem_phyid: np.ndarray       # (n_elements,) int, first tag
    elem_nodes: list             # list of 0-based node-index tuples
    physical_names: dict         # {(dim, phyid): name}

    @property
    def n_nodes(self) -> int:
        return self.nodes.shape[0]

    @property
    def n_elements(self) -> int:
        return self.elem_type.size


def read_gmsh22(path) -> GmshMesh:
    path = Path(path)
    lines = path.read_text().splitlines()
    i = 0
    n = len(lines)

    nodes = None
    node_ids = None
    elem_type: list[int] = []
    elem_phyid: list[int] = []
    elem_nodes: list[tuple] = []
    physical_names: dict = {}

    while i < n:
        tag = lines[i].strip()
        if tag == "$MeshFormat":
            ver = lines[i + 1].split()
            if not ver[0].startswith("2.2"):
                raise ValueError(f"{path}: expected gmsh ASCII 2.2, got {ver[0]}")
            if ver[1] != "0":
                raise ValueError(f"{path}: binary gmsh files are not supported")
            i += 2
        elif tag == "$PhysicalNames":
            cnt = int(lines[i + 1])
            for k in range(cnt):
                parts = lines[i + 2 + k].split(None, 2)
                physical_names[(int(parts[0]), int(parts[1]))] = parts[2].strip().strip('"')
            i += 2 + cnt
        elif tag == "$Nodes":
            cnt = int(lines[i + 1])
            nodes = np.zeros((cnt, 2), dtype=np.float64)
            node_ids = np.zeros(cnt, dtype=np.int64)
            for k in range(cnt):
                p = lines[i + 2 + k].split()
                node_ids[k] = int(p[0])
                nodes[k, 0] = float(p[1])
                nodes[k, 1] = float(p[2])
            i += 2 + cnt
        elif tag == "$Elements":
            cnt = int(lines[i + 1])
            for k in range(cnt):
                p = lines[i + 2 + k].split()
                etype = int(p[1])
                ntags = int(p[2])
                if etype not in (LINE_ELEMENT, TRIANGLE_ELEMENT):
                    raise ValueError(
                        f"{path}: element {p[0]} has unsupported type {etype}. "
                        "The Fortran reference silently desynchronises its file "
                        "pointer on such elements; refusing to guess.")
                if ntags < 2:
                    raise ValueError(f"{path}: element {p[0]} has {ntags} tags, need >= 2")
                nn = 2 if etype == LINE_ELEMENT else 3
                elem_type.append(etype)
                elem_phyid.append(int(p[3]))
                elem_nodes.append(tuple(int(v) for v in p[3 + ntags:3 + ntags + nn]))
            i += 2 + cnt
        else:
            i += 1

    if nodes is None:
        raise ValueError(f"{path}: no $Nodes section")

    # The Fortran indexes NODES by the *position* in the file, ignoring the
    # printed id.  Match that, but check the ids are 1..N so nothing silently
    # shifts on a renumbered mesh.
    expected = np.arange(1, nodes.shape[0] + 1)
    if not np.array_equal(node_ids, expected):
        raise ValueError(f"{path}: node ids are not 1..{nodes.shape[0]} in order; "
                         "the reference reader assumes they are")

    remap = {gid: k for k, gid in enumerate(node_ids)}
    elem_nodes0 = [tuple(remap[g] for g in e) for e in elem_nodes]

    return GmshMesh(nodes=nodes,
                    elem_type=np.asarray(elem_type, dtype=np.int64),
                    elem_phyid=np.asarray(elem_phyid, dtype=np.int64),
                    elem_nodes=elem_nodes0,
                    physical_names=physical_names)
