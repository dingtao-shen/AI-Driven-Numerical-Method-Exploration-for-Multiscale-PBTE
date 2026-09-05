"""Spatial mesh: gmsh reading, face topology, element geometry."""
from .geometry import triangle_area, triangle_geometry
from .gmsh_reader import GmshMesh, read_gmsh22
from .topology import Mesh, build_mesh

__all__ = ["GmshMesh", "read_gmsh22", "Mesh", "build_mesh",
           "triangle_area", "triangle_geometry"]
