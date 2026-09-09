"""pybte -- steady 2D2V gray linear Callaway phonon BTE, nodal DG + discrete ordinates.

A fidelity port of the Fortran research solver in ``fortran-reference/``,
offering the same two iteration schemes:

``accflag = 0``  CIS  -- conventional source iteration
``accflag = 1``  GSIS -- general synthetic iterative scheme (HDG macroscopic
                 acceleration)

Typical use::

    from pybte import Case, Solver

    case = Case.from_yaml("cases/cavity_tauR1e-3_gsis.yaml")
    record = Solver(case).run()
    record.iterations, record.converged, record.temp
"""
from .config import (Boundary, Case, DG, Flow, Iteration, MeshCfg, Output,
                     Performance, Restart, Scheme, VelMesh,
                     case_from_control_in, parse_control_in)
from .driver import Solver
from .io_output import RunRecord

__all__ = ["Case", "Solver", "RunRecord", "Boundary", "Iteration", "Scheme",
           "VelMesh", "DG", "Flow", "MeshCfg", "Restart", "Output",
           "Performance", "case_from_control_in", "parse_control_in"]

__version__ = "1.0.0"
