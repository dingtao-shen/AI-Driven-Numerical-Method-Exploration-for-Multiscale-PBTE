"""Per-direction sweep ordering -- the port of ``Matrix.f90::Init_Triangle_Order``.

For each discrete direction the elements are topologically sorted so that
every element's upwind neighbours are visited before it.  That is what turns
``DG_Solver_VDF`` into a genuine transport sweep (Gauss-Seidel in the sweep
order) rather than a Jacobi update, and it is why the ordering has to be
reproduced *exactly*, not merely validly:

    the Fortran repeatedly scans elements in ascending index order and emits
    every element whose incoming-face flags have all been cleared, clearing
    the reciprocal flag on its neighbours *immediately*, so an element later
    in the same scan can become ready and be emitted in the same pass.

Any other valid topological order would still converge, but to a different
iterate sequence, and every iteration count downstream would drift.

The Fortran's ``DO`` loop has no cycle detection and spins
forever if a direction produces a cyclic dependency.  Here a pass that emits
nothing raises :class:`SweepCycleError` by default; ``on_cycle='break'``
instead emits the element with the fewest incoming faces (lagging its inflow
to the previous iterate) and records which directions needed it.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._numba import njit, prange

__all__ = ["SweepOrder", "build_sweep_order", "SweepCycleError",
           "validate_sweep_order"]


class SweepCycleError(RuntimeError):
    """A direction's element dependency graph contains a cycle."""


@dataclass
class SweepOrder:
    """``order[d, k]`` is the k-th element to solve for direction ``d``.

    Directions are flattened in Fortran order over ``(J1, J2)``, i.e.
    ``d = j1 + npole * j2``, matching ``TRI_ORDER(:, J1, J2)``.
    """
    order: np.ndarray            # (ndir, n_tris) int32
    npole: int
    nazim: int
    cycles_broken: np.ndarray    # (ndir,) int32, elements force-emitted

    def fortran_tri_order(self) -> np.ndarray:
        """``TRI_ORDER(N_TRIS, NPOLE, NAZIM)`` with 1-based element ids."""
        n_tris = self.order.shape[1]
        return (self.order.reshape(self.nazim, self.npole, n_tris)
                .transpose(2, 1, 0) + 1).astype(np.int32)

    def __getitem__(self, d):
        return self.order[d]


@njit(cache=True)
def _order_kernel(neighbour, nx, ny, cx, cy, order, cycles, allow_break):
    """Topological sort per direction; see the module docstring.

    ``neighbour[i, il]`` is the element across local edge ``il`` or -1.
    ``nx``/``ny`` are ``(n_tris, 3)`` outward unit normals.
    ``cx``/``cy`` are ``(ndir,)``.  Returns 0 on success, -1 on a cycle.
    """
    n_tris = neighbour.shape[0]
    ndir = cx.shape[0]
    status = 0
    for d in prange(ndir):
        ioflag = np.zeros((n_tris, 3), dtype=np.int8)
        inarray = np.ones(n_tris, dtype=np.bool_)
        for i in range(n_tris):
            for il in range(3):
                if neighbour[i, il] >= 0:
                    flux = cx[d] * nx[i, il] + cy[d] * ny[i, il]
                    if flux > 0.0:
                        ioflag[i, il] = 1
                    elif flux < 0.0:
                        ioflag[i, il] = -1
        num = 0
        broken = 0
        while num < n_tris:
            progressed = False
            for i in range(n_tris):
                if not inarray[i]:
                    continue
                outgoing = True
                for il in range(3):
                    if ioflag[i, il] == -1:
                        outgoing = False
                        break
                if outgoing:
                    order[d, num] = i
                    num += 1
                    inarray[i] = False
                    progressed = True
                    for il in range(3):
                        t2 = neighbour[i, il]
                        if t2 >= 0:
                            for k in range(3):
                                if neighbour[t2, k] == i:
                                    ioflag[t2, k] = 0
                                    break
            if not progressed:
                if not allow_break:
                    status = -1
                    break
                # break the cycle: emit the still-queued element with the
                # fewest remaining incoming faces, lowest index wins
                best = -1
                best_in = 4
                for i in range(n_tris):
                    if not inarray[i]:
                        continue
                    cnt = 0
                    for il in range(3):
                        if ioflag[i, il] == -1:
                            cnt += 1
                    if cnt < best_in:
                        best_in = cnt
                        best = i
                if best < 0:
                    status = -1
                    break
                order[d, num] = best
                num += 1
                inarray[best] = False
                broken += 1
                for il in range(3):
                    t2 = neighbour[best, il]
                    if t2 >= 0:
                        for k in range(3):
                            if neighbour[t2, k] == best:
                                ioflag[t2, k] = 0
                                break
        cycles[d] = broken
    return status


def build_sweep_order(mesh, vel, on_cycle: str = "raise") -> SweepOrder:
    """Build ``TRI_ORDER`` for every direction.

    Parameters
    ----------
    mesh : :class:`pybte.mesh.Mesh`
    vel : :class:`pybte.velocity.VelocityMesh`
    on_cycle : {'raise', 'break'}
    """
    if on_cycle not in ("raise", "break"):
        raise ValueError("on_cycle must be 'raise' or 'break'")
    cx, cy, _ = vel.flat()
    ndir = cx.size
    n_tris = mesh.n_tris
    order = np.full((ndir, n_tris), -1, dtype=np.int32)
    cycles = np.zeros(ndir, dtype=np.int32)

    nx = np.ascontiguousarray(mesh.tri_normal[:, :, 0])
    ny = np.ascontiguousarray(mesh.tri_normal[:, :, 1])
    neighbour = np.ascontiguousarray(mesh.tri_neighbour.astype(np.int32))

    status = _order_kernel(neighbour, nx, ny,
                           np.ascontiguousarray(cx), np.ascontiguousarray(cy),
                           order, cycles, on_cycle == "break")
    if status != 0:
        bad = int(np.flatnonzero((order < 0).any(axis=1))[0])
        raise SweepCycleError(
            f"direction {bad} (j1={bad % vel.npole}, j2={bad // vel.npole}) has a "
            "cyclic element dependency; the mesh is not sweepable for it. "
            "Pass scheme.on_cycle='break' to lag the inflow instead.")

    return SweepOrder(order=order, npole=vel.npole, nazim=vel.nazim,
                      cycles_broken=cycles)


def validate_sweep_order(mesh, vel, sweep: SweepOrder) -> None:
    """Assert that every element's upwind neighbours precede it.

    This is the check the Fortran runs after building ``TRI_ORDER`` (the
    ``'Order Err'`` block), reproduced as a proper exception.
    """
    cx, cy, _ = vel.flat()
    for d in range(cx.size):
        pos = np.empty(mesh.n_tris, dtype=np.int64)
        pos[sweep.order[d]] = np.arange(mesh.n_tris)
        if sweep.cycles_broken[d]:
            continue
        for il in range(3):
            nb = mesh.tri_neighbour[:, il]
            has = nb >= 0
            flux = cx[d] * mesh.tri_normal[:, il, 0] + cy[d] * mesh.tri_normal[:, il, 1]
            idx = np.flatnonzero(has)
            downwind = idx[flux[idx] > 0]
            upwind = idx[flux[idx] < 0]
            if np.any(pos[nb[downwind]] < pos[downwind]):
                raise SweepCycleError(f"outflow neighbour precedes element, direction {d}")
            if np.any(pos[nb[upwind]] > pos[upwind]):
                raise SweepCycleError(f"inflow neighbour follows element, direction {d}")
