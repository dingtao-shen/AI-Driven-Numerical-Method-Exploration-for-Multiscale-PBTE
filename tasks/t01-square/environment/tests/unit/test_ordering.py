"""The per-direction sweep ordering, cycle handling included."""
from __future__ import annotations

import numpy as np
import pytest

from pybte.ordering import (SweepCycleError, build_sweep_order,
                            validate_sweep_order)
from pybte.velocity import build_velocity_mesh


@pytest.fixture(scope="module")
def order(shipped_mesh):
    return build_sweep_order(shipped_mesh, build_velocity_mesh(20, 40, 1.0))


def test_order_is_a_permutation(order, shipped_mesh):
    n = shipped_mesh.n_tris
    for d in range(order.order.shape[0]):
        assert np.array_equal(np.sort(order.order[d]), np.arange(n))


def test_order_is_topologically_valid(shipped_mesh, order):
    """Every element's upwind neighbours precede it -- this is what makes the
    sweep a true Gauss-Seidel transport solve rather than a Jacobi update."""
    validate_sweep_order(shipped_mesh, build_velocity_mesh(20, 40, 1.0), order)


def test_no_cycles_on_the_shipped_mesh(order):
    assert int(order.cycles_broken.sum()) == 0


def test_tie_breaking_is_ascending_index(shipped_mesh, order):
    """The Fortran scans elements in ascending index order and emits every
    ready one immediately, clearing neighbour flags as it goes.  Reproduce
    that here with a plain Python implementation and require equality: any
    other valid topological order would converge differently."""
    vel = build_velocity_mesh(20, 40, 1.0)
    cx, cy, _ = vel.flat()
    m = shipped_mesh
    n = m.n_tris
    for d in (0, 137, 799):
        flag = np.zeros((n, 3), dtype=int)
        for i in range(n):
            for il in range(3):
                if m.tri_neighbour[i, il] >= 0:
                    f = cx[d] * m.tri_normal[i, il, 0] + cy[d] * m.tri_normal[i, il, 1]
                    flag[i, il] = 1 if f > 0 else (-1 if f < 0 else 0)
        inarray = np.ones(n, dtype=bool)
        out = []
        while len(out) < n:
            progressed = False
            for i in range(n):
                if not inarray[i]:
                    continue
                if np.any(flag[i] == -1):
                    continue
                out.append(i)
                inarray[i] = False
                progressed = True
                for il in range(3):
                    t2 = m.tri_neighbour[i, il]
                    if t2 >= 0:
                        for k in range(3):
                            if m.tri_neighbour[t2, k] == i:
                                flag[t2, k] = 0
                                break
            assert progressed
        assert np.array_equal(np.asarray(out), order.order[d])


def test_cycle_detection_raises(monkeypatch, shipped_mesh):
    """The Fortran's loop spins forever on a cyclic direction.  We
    construct one by lying about the neighbour graph and check we raise."""
    import copy

    m = copy.copy(shipped_mesh)
    nb = shipped_mesh.tri_neighbour.copy()
    # make a 3-cycle: force elements 0,1,2 to point at each other as inflow
    # regardless of geometry by giving them identical inward normals
    normal = shipped_mesh.tri_normal.copy()
    nb[0] = [1, -1, -1]
    nb[1] = [2, -1, -1]
    nb[2] = [0, -1, -1]
    normal[0, 0] = [-1.0, 0.0]
    normal[1, 0] = [-1.0, 0.0]
    normal[2, 0] = [-1.0, 0.0]
    m.tri_neighbour = nb
    m.tri_normal = normal
    vel = build_velocity_mesh(4, 4, 1.0)
    with pytest.raises(SweepCycleError):
        build_sweep_order(m, vel, on_cycle="raise")


def test_cycle_breaking_completes(shipped_mesh):
    import copy

    m = copy.copy(shipped_mesh)
    nb = shipped_mesh.tri_neighbour.copy()
    normal = shipped_mesh.tri_normal.copy()
    nb[0] = [1, -1, -1]
    nb[1] = [2, -1, -1]
    nb[2] = [0, -1, -1]
    normal[0, 0] = [-1.0, 0.0]
    normal[1, 0] = [-1.0, 0.0]
    normal[2, 0] = [-1.0, 0.0]
    m.tri_neighbour = nb
    m.tri_normal = normal
    so = build_sweep_order(m, build_velocity_mesh(4, 4, 1.0), on_cycle="break")
    assert int(so.cycles_broken.sum()) > 0
    for d in range(so.order.shape[0]):
        assert np.array_equal(np.sort(so.order[d]), np.arange(m.n_tris))
