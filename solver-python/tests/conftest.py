"""Shared fixtures for the pybte test suite."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

MESH = ROOT / "meshes" / "A1_Nx11_Ny11.msh"
CASES = ROOT / "cases"


@pytest.fixture(scope="session")
def shipped_case():
    from pybte import Case
    return Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")


@pytest.fixture(scope="session")
def shipped_mesh():
    from pybte.mesh import build_mesh
    return build_mesh(MESH, [11, 12, 13, 14], [1, 1, 1, 1],
                      ["SWall", "NWall", "EWall", "WWall"])


@pytest.fixture(scope="session")
def basis3():
    from pybte.basis import build_basis
    return build_basis(3)
