"""Paths shared by the unit tests that need to build a Case by hand."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASES = ROOT / "cases"
MESH_COARSE = ROOT / "meshes" / "A1_Nx6_Ny6.msh"
