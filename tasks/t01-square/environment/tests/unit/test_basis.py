"""Nodal bases."""
from __future__ import annotations

import numpy as np
import pytest

from pybte.basis import build_basis, eval_tri_basis, monomial_index


@pytest.mark.parametrize("deg", [1, 2, 3, 4])
def test_tri_basis_kronecker_delta(deg):
    b = build_basis(deg)
    v = eval_tri_basis(b.nodfun_tri, b.mono, b.node_tri[:, 0], b.node_tri[:, 1])
    assert np.allclose(v, np.eye(b.ndof_tri), atol=1e-12)


@pytest.mark.parametrize("deg", [1, 2, 3, 4])
def test_tri_basis_partition_of_unity(deg):
    b = build_basis(deg)
    rng = np.random.default_rng(1234)
    xi = rng.random(400)
    eta = rng.random(400) * (1.0 - xi)
    v = eval_tri_basis(b.nodfun_tri, b.mono, xi, eta)
    assert np.allclose(v.sum(axis=1), 1.0, atol=1e-11)


@pytest.mark.parametrize("deg", [1, 2, 3, 4])
def test_tri_basis_reproduces_polynomials(deg):
    """Interpolating a polynomial of degree <= deg is exact."""
    b = build_basis(deg)
    rng = np.random.default_rng(7)
    coeffs = rng.normal(size=b.ndof_tri)

    def f(x, y):
        return sum(c * x ** float(a) * y ** float(bb)
                   for c, (a, bb) in zip(coeffs, b.mono))

    nodal = f(b.node_tri[:, 0], b.node_tri[:, 1])
    xi = rng.random(200)
    eta = rng.random(200) * (1.0 - xi)
    got = eval_tri_basis(b.nodfun_tri, b.mono, xi, eta) @ nodal
    assert np.allclose(got, f(xi, eta), atol=1e-10)


@pytest.mark.parametrize("deg", [1, 2, 3, 4])
def test_fc_basis(deg):
    b = build_basis(deg)
    t = b.node_fc
    v = b.eval_fc(t)
    assert np.allclose(v, np.eye(b.ndof_fc), atol=1e-12)
    tt = np.linspace(-1, 1, 37)
    assert np.allclose(b.eval_fc(tt).sum(axis=1), 1.0, atol=1e-11)


def test_monomial_index_matches_ordering():
    b = build_basis(3)
    for k, (a, c) in enumerate(b.mono):
        assert monomial_index(int(a), int(c)) == k


def test_node_layout():
    """Nodes run j = 0..DEG outer, i = 0..DEG-j inner, at (i/DEG, j/DEG)."""
    b = build_basis(3)
    expect = [(i / 3, j / 3) for j in range(4) for i in range(4 - j)]
    assert np.allclose(b.node_tri, np.array(expect))
