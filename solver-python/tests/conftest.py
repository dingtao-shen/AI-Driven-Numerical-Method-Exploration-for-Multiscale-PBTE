"""Shared fixtures and skip logic for the pybte test suite."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

GOLDEN = REPO / "fortran-reference" / "golden"
MESH = ROOT / "meshes" / "A1_Nx11_Ny11.msh"
CASES = ROOT / "cases"


def has_dump(name: str) -> bool:
    return (GOLDEN / name / "dump" / "manifest.txt").exists()


def has_golden(name: str) -> bool:
    return (GOLDEN / name / "residual_history.txt").exists()


requires_dump = pytest.mark.skipif(
    not has_dump("dump_cis_shipped"),
    reason="Fortran stage dumps not built; run tools/build_fortran.sh and "
           "tools/make_golden.py first")


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
