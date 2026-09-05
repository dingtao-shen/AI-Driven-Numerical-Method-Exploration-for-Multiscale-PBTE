#!/usr/bin/env python3
"""Generate the gmsh 2.2 ASCII meshes used by the cases and the test suite.

The shipped ``A1_Nx11_Ny11.msh`` is a structured triangulation of the unit
square whose node coordinates follow the smootherstep clustering

    s(t) = t^3 (10 - 15 t + 6 t^2),      t = i/N

which packs points into the two boundary layers -- exactly what a kinetic
problem at small Knudsen number needs.  Reverse-engineering it lets us build
a refinement *family* with the same character, so mesh convergence studies
compare like with like.  ``--check`` verifies that the N=10 member reproduces
the shipped file byte for byte.

Layouts
-------
``cavity``   four walls, physical ids 11 (south), 12 (north), 13 (east),
             14 (west) -- the shipped configuration
``channel``  north/south walls (11, 12) plus a periodic pair on east/west
             named ``Master`` (13, x=1) and ``Slave`` (14, x=0)

Usage::

    python tools/make_meshes.py --outdir meshes
    python tools/make_meshes.py --check meshes/A1_Nx11_Ny11.msh
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def smootherstep(t):
    t = np.asarray(t, dtype=float)
    return t ** 3 * (10.0 - 15.0 * t + 6.0 * t ** 2)


def uniform(t):
    return np.asarray(t, dtype=float)


DISTRIBUTIONS = {"smootherstep": smootherstep, "uniform": uniform}

_HEADER = """$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
5
1 11 "{s} (physical id 11)"
1 12 "{n} (physical id 12)"
1 13 "{e} (physical id 13)"
1 14 "{w} (physical id 14)"
2 15 "Fluid (physical id 15)"
$EndPhysicalNames
"""


def build(n: int, layout: str = "cavity", distribution: str = "smootherstep",
          distribution_x: str | None = None) -> str:
    """Return the file contents for an ``n x n`` cell grid (``(n+1)^2`` nodes).

    ``distribution_x`` defaults to ``distribution``.  The channel uses
    ``uniform`` in x so that the mesh is genuinely invariant under the
    periodic translation: a graded x-spacing makes neighbouring columns
    non-congruent, and the converged field then carries a few-percent
    x-dependence that is a property of the mesh, not of the boundary
    condition.
    """
    fy = DISTRIBUTIONS[distribution]
    fx = DISTRIBUTIONS[distribution_x or distribution]
    coord_y = fy(np.arange(n + 1) / n)
    coord_x = fx(np.arange(n + 1) / n)
    nn = n + 1

    names = {"s": "Swall", "n": "NWall", "e": "EWall", "w": "WWall"}
    if layout == "channel":
        names = {"s": "Swall", "n": "NWall", "e": "Master", "w": "Slave"}
    elif layout != "cavity":
        raise ValueError(f"unknown layout {layout!r}")

    out = [_HEADER.format(**names)]

    def nid(i, j):
        """1-based node id; the shipped file numbers x-major."""
        return i * nn + j + 1

    out.append("$Nodes\n")
    out.append(f"{nn * nn:6d}\n")
    for i in range(nn):
        for j in range(nn):
            out.append(f"{nid(i, j):6d} {coord_x[i]:.16f} {coord_y[j]:.16f} 0\n")
    out.append("$EndNodes\n")

    lines = []
    for i in range(n):                       # south, +x along y=0
        lines.append((11, 5, nid(i, 0), nid(i + 1, 0)))
    for i in range(n, 0, -1):                # north, -x along y=1
        lines.append((12, 6, nid(i, n), nid(i - 1, n)))
    for j in range(n):                       # east, +y along x=1
        lines.append((13, 7, nid(n, j), nid(n, j + 1)))
    for j in range(n, 0, -1):                # west, -y along x=0
        lines.append((14, 8, nid(0, j), nid(0, j - 1)))

    # The shipped mesh chooses the quad diagonal per quadrant: it flips
    # across x = 1/2 and again across y = 1/2, so the diagonal always points
    # towards the domain centre -- the shorter one for this graded spacing,
    # and symmetric about both mid-planes.  Reproduce it exactly: the
    # triangulation fixes the element ordering, and the element ordering is
    # TRI_ORDER's tie-breaker.
    # ... except for the channel, where the east and west boundaries are
    # identified.  A per-quadrant flip makes the triangulation on the two
    # sides of the periodic interface *different*, so the discrete problem is
    # not translationally invariant and the converged field acquires a
    # spurious few-percent x-dependence.  For the channel the diagonal
    # therefore flips on y only, and every x-column is identical.
    flip_on_x = layout != "channel"

    tris = []
    for i in range(n):
        for j in range(n):
            if (2 * i >= n if flip_on_x else False) == (2 * j >= n):
                tris.append((nid(i, j), nid(i + 1, j), nid(i + 1, j + 1)))
                tris.append((nid(i, j), nid(i + 1, j + 1), nid(i, j + 1)))
            else:
                tris.append((nid(i, j), nid(i + 1, j), nid(i, j + 1)))
                tris.append((nid(i, j + 1), nid(i + 1, j), nid(i + 1, j + 1)))

    out.append("$Elements\n")
    out.append(f"{len(lines) + len(tris):6d}\n")
    k = 0
    for phy, geo, a, b in lines:
        k += 1
        out.append(f"{k:6d} 1 2 {phy} {geo} {a:5d} {b:5d}\n")
    for a, b, c in tris:
        k += 1
        out.append(f"{k:6d} 2 2 15 10 {a:5d} {b:5d} {c:5d}\n")
    out.append("$EndElements\n")
    return "".join(out)


def _numbers(text):
    """Every numeric token, for a format-insensitive comparison."""
    toks = []
    for line in text.splitlines():
        if line.startswith("$") or '"' in line:
            continue
        for t in line.split():
            try:
                toks.append(float(t))
            except ValueError:
                pass
    return np.asarray(toks)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="meshes")
    p.add_argument("--levels", default="5,10,20,40",
                   help="cells per side for the refinement family")
    p.add_argument("--check", default=None,
                   help="verify the generator reproduces this shipped mesh")
    args = p.parse_args(argv)

    if args.check:
        shipped = Path(args.check).read_text()
        mine = build(10, "cavity", "smootherstep")
        a, b = _numbers(shipped), _numbers(mine)
        if a.shape != b.shape:
            print(f"FAIL: {a.shape} vs {b.shape} numeric tokens")
            return 1
        d = np.abs(a - b).max()
        print(f"{'OK' if d < 1e-12 else 'FAIL'}: max token difference {d:.3e}")
        return 0 if d < 1e-12 else 1

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    for n in (int(v) for v in args.levels.split(",")):
        for layout, tag, dx in (("cavity", "A", None), ("channel", "C", "uniform")):
            path = out / f"{tag}1_Nx{n + 1}_Ny{n + 1}.msh"
            path.write_text(build(n, layout, distribution_x=dx))
            print(f"wrote {path}  ({2 * n * n} triangles)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
