"""Constants mirroring ``Global_Parameter.f90::MODULE CONSTANT``.

The Fortran uses ``PI = 4*atan(1)`` in double precision, which is bit-for-bit
``math.pi``.  ``DBL`` is IEEE binary64, i.e. numpy ``float64``.
"""
from __future__ import annotations

import numpy as np

#: ``REAL(KIND=DBL)`` -- the only floating type used anywhere in the solver.
DBL = np.float64

PI: float = 4.0 * np.arctan(1.0)

BOLTZ: float = 1.3806485279e-23
AFGDL: float = 6.022140857e23

MIN_ERR: float = 1.0e-16
MID_ERR: float = 1.0e-10
MAX_ERR: float = 1.0e-6

#: Boundary-condition type codes, as used by ``BC_TYP`` in ``control.in``.
BC_THERMALISING = 1      # isothermal / diffusely emitting wall
BC_NONTHERMALISING = 2   # adiabatic, diffusely reflecting wall
BC_PERIODIC = 3          # periodic pair
BC_SYMMETRY = 4          # specular symmetry plane (mesh bookkeeping only)

BC_TYPE_NAMES = {
    "thermalising": BC_THERMALISING,
    "thermalizing": BC_THERMALISING,
    "isothermal": BC_THERMALISING,
    "nonthermalising": BC_NONTHERMALISING,
    "nonthermalizing": BC_NONTHERMALISING,
    "adiabatic": BC_NONTHERMALISING,
    "periodic": BC_PERIODIC,
    "symmetry": BC_SYMMETRY,
}
BC_NAME_BY_CODE = {
    BC_THERMALISING: "thermalising",
    BC_NONTHERMALISING: "nonthermalising",
    BC_PERIODIC: "periodic",
    BC_SYMMETRY: "symmetry",
}

#: ``NP_TRI`` as a function of ``DEG`` (``Init_Global_Variables``).
NP_TRI_BY_DEG = {1: 3, 2: 6, 3: 12, 4: 16}

#: ``NP_FC`` is hard-coded to 15 in the Fortran.
NP_FC_DEFAULT = 15
