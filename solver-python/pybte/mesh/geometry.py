"""Element geometry -- signed area, outward edge normals, inscribed heights.

Transcribed from the second half of ``Spatial_Mesh.f90::Init_Spatial_Grid``.
The floating-point expressions are written in the Fortran's exact order so
that the results are bit-identical, not merely equal to rounding.

Edge convention (``Spatial_Mesh.f90`` header comment)::

      3
      | \\        edge 1: node 1 -- node 2
      |  \\       edge 2: node 2 -- node 3
      1---2      edge 3: node 3 -- node 1

For a counter-clockwise triangle the outward normal of edge ``j`` is
``(dy, -dx)/|d|`` with ``d`` the edge vector, which is what the Fortran
writes as ``((y_{j+1}-y_j)/L, (x_j-x_{j+1})/L)``.
"""
from __future__ import annotations

import numpy as np

__all__ = ["triangle_area", "triangle_geometry"]


def triangle_area(x, y) -> np.ndarray:
    """Signed area, in the Fortran's exact summation order.

    ``0.5*(x3*y1 - x2*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3)``
    """
    x1, x2, x3 = x[:, 0], x[:, 1], x[:, 2]
    y1, y2, y3 = y[:, 0], y[:, 1], y[:, 2]
    return 0.5 * (x3 * y1 - x2 * y1 + x1 * y2 - x3 * y2 - x1 * y3 + x2 * y3)


def triangle_geometry(nodes: np.ndarray, tri_nodes: np.ndarray):
    """Return ``(area, normal, hmin_per_tri, hmin_global)``.

    Parameters
    ----------
    nodes : (n_nodes, 2)
    tri_nodes : (n_tris, 3) int, 0-based

    Returns
    -------
    area : (n_tris,)
    normal : (n_tris, 3, 2) -- ``normal[i, j]`` is the outward unit normal of
        local edge ``j`` of triangle ``i``
    hmin_per_tri : (n_tris,) -- the smallest of the three heights ``2A/L_j``
    hmin_global : float -- ``min(1.0, min_ij 2A_i/L_ij)``, exactly as the
        Fortran initialises ``Hmin = 1.d0`` before the scan
    """
    n_tris = tri_nodes.shape[0]
    x = nodes[tri_nodes, 0]                       # (n_tris, 3)
    y = nodes[tri_nodes, 1]

    area = triangle_area(x, y)
    if np.any(area < 0):
        bad = np.flatnonzero(area < 0)
        raise ValueError(f"{bad.size} triangles have negative area (node ordering); "
                         f"first offenders: {bad[:5].tolist()}")

    xc = np.concatenate([x, x[:, :1]], axis=1)    # x(4) = x(1)
    yc = np.concatenate([y, y[:, :1]], axis=1)

    dx = xc[:, :3] - xc[:, 1:4]                   # nx(J) - nx(J+1)
    dy = yc[:, :3] - yc[:, 1:4]
    length = np.sqrt(dx * dx + dy * dy)

    normal = np.empty((n_tris, 3, 2), dtype=np.float64)
    normal[:, :, 0] = -dy / length                # (ny(J+1)-ny(J))/L
    normal[:, :, 1] = dx / length                 # (nx(J)-nx(J+1))/L

    # sl uses **2 in the Fortran rather than the product form; same value
    sl = np.sqrt(dx ** 2 + dy ** 2)
    hl = 2.0 * area[:, None] / sl
    hmin_per_tri = hl.min(axis=1)
    hmin_global = min(1.0, float(hl.min()))

    return area, normal, hmin_per_tri, hmin_global
