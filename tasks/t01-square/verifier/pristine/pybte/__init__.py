"""pybte -- steady 2D2V gray linear Callaway phonon BTE, nodal DG + discrete ordinates.

Typical use::

    from pybte import Case, Solver

    case = Case.from_yaml("cases/kn_1e-2.yaml")
    record = Solver(case).run()
    record.iterations, record.converged, record.temp
"""
from .config import (Boundary, Case, DG, Flow, Iteration, MeshCfg, Output,
                     Performance, Restart, Scheme, VelMesh)
from .driver import Solver
from .io_output import RunRecord

__all__ = ["Case", "Solver", "RunRecord", "Boundary", "Iteration", "Scheme",
           "VelMesh", "DG", "Flow", "MeshCfg", "Restart", "Output",
           "Performance"]

__version__ = "1.0.0"
