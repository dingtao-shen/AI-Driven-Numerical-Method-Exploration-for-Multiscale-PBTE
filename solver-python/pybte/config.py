"""Case configuration -- the YAML replacement for ``control.in``.

The mapping is one-to-one with the Fortran namelists,
plus the fields that were implicit in the Fortran and are made explicit here:

* ``scheme.acc_variant``   -- A (``Synthetic_Acceleration.f90``, the default)
                              or B (``Synthetic_Acceleration1.f90``)
* ``scheme.stabilisation`` -- the HDG ``ST(1:3)``
* ``restart``              -- explicit, off by default
* ``iteration.true_residual`` -- report a true transport residual too
* ``scheme.method``        -- ``source`` (default) or ``krylov``: GMRES on
                              the outer iteration, same fixed point as CIS;
                              see docs/KRYLOV.md
* ``scheme.defect_*``      -- experimental fixed-point repair, off by
                              default; see docs/DEFECT_CORRECTION.md
* ``performance``          -- storage/kernel choices, no physics
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any

from .constants import (BC_NAME_BY_CODE, BC_TYPE_NAMES, NP_FC_DEFAULT,
                        NP_TRI_BY_DEG)

__all__ = ["Case", "Boundary", "Iteration", "Scheme", "VelMesh", "DG", "Flow",
           "MeshCfg", "Restart", "Output", "Performance", "parse_control_in"]


# --------------------------------------------------------------------------
@dataclass
class Iteration:
    tol: float = 1.0e-8
    tmax: int = 10_000
    #: Also compute the transport residual of the discrete system each
    #: iteration.  Never changes the stopping criterion, which stays
    #: ``residual_iterate < tol``.  For CIS this tracks the iterate residual;
    #: what exposes the error is ``RunRecord.error_estimate``.
    true_residual: bool = False
    #: evaluate the true residual only every k-th iteration (it costs a sweep)
    true_residual_every: int = 1

    def __post_init__(self):
        # YAML 1.1 reads "1.0e5" (no exponent sign) as a string, and the
        # Fortran namelist habit is to write exactly that.  Coerce rather
        # than fail three call levels later with a TypeError.
        self.tol = float(self.tol)
        self.tmax = int(self.tmax)
        self.true_residual = bool(self.true_residual)
        self.true_residual_every = int(self.true_residual_every)


@dataclass
class Scheme:
    accflag: int = 0                     # 0 = CIS, 1 = GSIS
    acc_variant: str = "A"
    stabilisation: list = field(default_factory=lambda: [1.0, 1.0, 1.0])
    #: divide ST(2:3) by the global minimum element height
    scale_by_hmin: bool = False
    #: 'raise' or 'break' when a direction's sweep graph has a cycle
    on_cycle: str = "raise"
    #: Experimental: relaxation on the lagged consistency defect that removes
    #: the GSIS/CIS fixed-point displacement.  ``0.0`` (default) is the
    #: shipped damped blend; see docs/DEFECT_CORRECTION.md.
    defect_omega: float = 0.0
    #: Update the defect only every N iterations.  Between updates it is
    #: frozen, which leaves the GSIS iteration matrix -- and therefore the
    #: convergence rate -- exactly as published; only the fixed point moves.
    defect_every: int = 1
    #: Anderson-acceleration window on the outer defect sequence.  0 = plain
    #: relaxation.  The outer map ``d -> o(d)`` is affine, so Anderson applies
    #: exactly; it contracts at ~0.95 unaccelerated, which is what makes the
    #: unaccelerated outer loop useless.
    defect_anderson: int = 0
    #: How the outer iteration is driven.  ``source``: Richardson / source
    #: iteration, the historical behaviour (with ``accflag = 1``, GSIS).
    #: ``krylov``: GMRES on ``(I - T) u = g``, one sweep per product -- same
    #: fixed point as source iteration; raises where a sweep is not an exact
    #: evaluation of the map (broken cycles, periodic faces).  See krylov.py.
    method: str = "source"
    krylov_tol_factor: float = 1.0e-6
    krylov_max_bytes: float = 5.0e8
    krylov_max_rounds: int = 3
    krylov_precond: bool = False

    def __post_init__(self):
        self.method = str(self.method).lower()
        if self.method not in ("source", "krylov"):
            raise ValueError("scheme.method must be 'source' or 'krylov'")
        self.krylov_tol_factor = float(self.krylov_tol_factor)
        self.krylov_max_bytes = float(self.krylov_max_bytes)
        self.krylov_max_rounds = int(self.krylov_max_rounds)
        self.krylov_precond = bool(self.krylov_precond)
        if self.krylov_max_bytes <= 0 or self.krylov_max_rounds < 1:
            raise ValueError("scheme.krylov_max_bytes must be > 0 and "
                             "krylov_max_rounds >= 1")
        self.accflag = int(self.accflag)
        self.acc_variant = str(self.acc_variant).upper()
        self.stabilisation = [float(v) for v in self.stabilisation]
        self.scale_by_hmin = bool(self.scale_by_hmin)
        self.defect_omega = float(self.defect_omega)
        self.defect_every = int(self.defect_every)
        if self.defect_every < 0:
            raise ValueError("scheme.defect_every must be >= 0")
        self.defect_anderson = int(self.defect_anderson)
        if self.defect_anderson < 0:
            raise ValueError("scheme.defect_anderson must be >= 0")
        if not 0.0 <= self.defect_omega <= 1.0:
            raise ValueError("scheme.defect_omega must be in [0, 1]")
        if self.on_cycle not in ("raise", "break"):
            raise ValueError("scheme.on_cycle must be 'raise' or 'break'")


@dataclass
class VelMesh:
    npole: int = 20
    nazim: int = 40

    def __post_init__(self):
        self.npole = int(self.npole)
        self.nazim = int(self.nazim)


@dataclass
class DG:
    deg: int = 3

    def __post_init__(self):
        self.deg = int(self.deg)


@dataclass
class Flow:
    cv: float = 1.0
    vg: float = 1.0
    tau_r: float = 1.0e-3
    tau_n: float = 1.0e5
    tau_thr: float = 1.0

    def __post_init__(self):
        self.cv = float(self.cv)
        self.vg = float(self.vg)
        self.tau_r = float(self.tau_r)
        self.tau_n = float(self.tau_n)
        self.tau_thr = float(self.tau_thr)
        if self.tau_r <= 0.0 or self.tau_n <= 0.0:
            raise ValueError("tau_r and tau_n must be positive")
        if self.vg <= 0.0:
            raise ValueError("vg must be positive")


@dataclass
class MeshCfg:
    file: str = "meshes/A1_Nx11_Ny11.msh"


@dataclass
class Boundary:
    name: str
    phyid: int
    type: str = "thermalising"
    temp: float = 0.0
    xoff: float = 0.0
    yoff: float = 0.0

    def __post_init__(self):
        self.phyid = int(self.phyid)
        self.temp = float(self.temp)
        self.xoff = float(self.xoff)
        self.yoff = float(self.yoff)

    @property
    def code(self) -> int:
        key = str(self.type).strip().lower()
        if key not in BC_TYPE_NAMES:
            raise ValueError(f"unknown boundary type {self.type!r}; "
                             f"expected one of {sorted(set(BC_TYPE_NAMES))}")
        return BC_TYPE_NAMES[key]


@dataclass
class Restart:
    enabled: bool = False
    path: str | None = None
    #: with restart disabled, refuse to run if a stale VDF*.out is lying about
    error_on_stale: bool = True


@dataclass
class Output:
    dir: str = "out"
    field: bool = True
    runtime_log: bool = True
    run_record: bool = True
    tecplot: bool = False
    npz: bool = True
    field_nx: int = 109
    field_ny: int = 109


@dataclass
class Performance:
    """Storage and kernel choices.  None of these changes a single bit of the
    answer -- they are checked by ``tests/unit/test_sweep_kernels.py``."""

    #: LU-factorise A_SOL once per (element, direction) instead of every
    #: iteration.  Costs N_TRIS*NDIR*NDOF^2*8 bytes (128 MB on the shipped
    #: case); turn it off for large meshes.
    precompute_inverse: bool = True
    #: "numba" -- the jitted sweep; "numpy" -- ``sweep_reference``, which is
    #: orders of magnitude slower and exists to pin down what the jitted
    #: kernels must reproduce.  Falls back to "numpy" if numba is absent.
    kernel: str = "numba"
    #: Threads for the jitted sweep.  Directions are independent, so any
    #: value gives identical results; the default of 0 means "let numba
    #: decide" (all cores).  Set it to 1 to leave the machine free, or when
    #: several solvers run side by side and would otherwise oversubscribe.
    threads: int = 0

    def __post_init__(self):
        self.precompute_inverse = bool(self.precompute_inverse)
        self.kernel = str(self.kernel).lower()
        if self.kernel not in ("numba", "numpy"):
            raise ValueError("performance.kernel must be 'numba' or 'numpy'")
        self.threads = int(self.threads)
        if self.threads < 0:
            raise ValueError("performance.threads must be >= 0")


# --------------------------------------------------------------------------
@dataclass
class Case:
    iteration: Iteration = field(default_factory=Iteration)
    scheme: Scheme = field(default_factory=Scheme)
    velmesh: VelMesh = field(default_factory=VelMesh)
    dg: DG = field(default_factory=DG)
    flow: Flow = field(default_factory=Flow)
    mesh: MeshCfg = field(default_factory=MeshCfg)
    boundaries: list = field(default_factory=list)
    restart: Restart = field(default_factory=Restart)
    output: Output = field(default_factory=Output)
    performance: Performance = field(default_factory=Performance)
    #: directory that relative paths in this case resolve against
    basedir: str = "."

    # -- derived quantities (Init_Global_Variables) ------------------------
    def __post_init__(self):
        # NAZIM = (NAZIM/2)*2 -- the Fortran forces it even
        self.velmesh.nazim = (self.velmesh.nazim // 2) * 2
        if self.velmesh.nazim < 2:
            raise ValueError("nazim must be >= 2 after the even-ing step")
        if self.dg.deg not in NP_TRI_BY_DEG:
            raise ValueError(f"deg must be one of {sorted(NP_TRI_BY_DEG)}")
        if self.scheme.acc_variant not in ("A", "B"):
            raise ValueError("scheme.acc_variant must be 'A' or 'B'")
        if len(self.scheme.stabilisation) != 3:
            raise ValueError("scheme.stabilisation must have three entries")

    @property
    def ndof_tri(self) -> int:
        return (self.dg.deg + 1) * (self.dg.deg + 2) // 2

    @property
    def ndof_fc(self) -> int:
        return self.dg.deg + 1

    @property
    def np_tri(self) -> int:
        return NP_TRI_BY_DEG[self.dg.deg]

    @property
    def np_fc(self) -> int:
        return NP_FC_DEFAULT

    @property
    def tau_c(self) -> float:
        return 1.0 / (1.0 / self.flow.tau_r + 1.0 / self.flow.tau_n)

    @property
    def nbc(self) -> int:
        return len(self.boundaries)

    def bc_codes(self) -> list:
        return [b.code for b in self.boundaries]

    def mesh_path(self) -> Path:
        p = Path(self.mesh.file)
        return p if p.is_absolute() else (Path(self.basedir) / p)

    def output_path(self) -> Path:
        p = Path(self.output.dir)
        return p if p.is_absolute() else (Path(self.basedir) / p)

    # -- (de)serialisation -------------------------------------------------
    def to_dict(self) -> dict:
        d = asdict(self)
        d.pop("basedir", None)
        return d

    def canonical_json(self) -> str:
        """Everything that can change the numbers, and nothing else."""
        d = self.to_dict()
        d.pop("output", None)
        d.pop("performance", None)
        return json.dumps(d, sort_keys=True, separators=(",", ":"))

    @property
    def config_hash(self) -> str:
        h = hashlib.sha256(self.canonical_json().encode()).hexdigest()
        mesh = self.mesh_path()
        if mesh.exists():
            h2 = hashlib.sha256(mesh.read_bytes()).hexdigest()
            h = hashlib.sha256((h + h2).encode()).hexdigest()
        return h

    def replace(self, **kw) -> "Case":
        return replace(self, **kw)

    # -- constructors ------------------------------------------------------
    @classmethod
    def from_dict(cls, data: dict, basedir: str = ".") -> "Case":
        def sub(klass, key):
            return klass(**(data.get(key) or {}))

        bnds = [Boundary(**b) for b in (data.get("boundaries") or [])]
        return cls(
            iteration=sub(Iteration, "iteration"),
            scheme=sub(Scheme, "scheme"),
            velmesh=sub(VelMesh, "velmesh"),
            dg=sub(DG, "dg"),
            flow=sub(Flow, "flow"),
            mesh=sub(MeshCfg, "mesh"),
            boundaries=bnds,
            restart=sub(Restart, "restart"),
            output=sub(Output, "output"),
            performance=sub(Performance, "performance"),
            basedir=basedir,
        )

    @classmethod
    def from_yaml(cls, path) -> "Case":
        import yaml

        path = Path(path)
        data = yaml.safe_load(path.read_text()) or {}
        return cls.from_dict(data, basedir=str(path.parent))

    def to_yaml(self, path) -> None:
        import yaml

        Path(path).write_text(yaml.safe_dump(self.to_dict(), sort_keys=False))


# --------------------------------------------------------------------------
# control.in -> Case
# --------------------------------------------------------------------------
_NML_RE = re.compile(r"&(\w+)(.*?)^\s*/\s*$", re.S | re.M)


def _fortran_scalar(tok: str):
    tok = tok.strip().rstrip(",").strip()
    if tok.startswith("'") or tok.startswith('"'):
        return tok[1:-1]
    t = tok.lower().replace("d", "e")
    try:
        if re.fullmatch(r"[+-]?\d+", t):
            return int(t)
        return float(t)
    except ValueError:
        return tok


def _split_values(rhs: str):
    """Split a namelist RHS on commas that are not inside quotes."""
    out, cur, quote = [], "", None
    for ch in rhs:
        if quote:
            cur += ch
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            cur += ch
        elif ch == ",":
            out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return [v for v in (s.strip() for s in out) if v]


def parse_control_in(path) -> dict:
    """Parse the Fortran namelist file into ``{GROUP: {KEY: value_or_list}}``."""
    text = Path(path).read_text()
    # strip trailing "!" comments outside quotes
    cleaned = []
    for line in text.splitlines():
        out, quote = "", None
        for ch in line:
            if quote:
                out += ch
                if ch == quote:
                    quote = None
            elif ch in "'\"":
                quote = ch
                out += ch
            elif ch == "!":
                break
            else:
                out += ch
        cleaned.append(out)
    text = "\n".join(cleaned)

    groups: dict[str, dict[str, Any]] = {}
    for m in _NML_RE.finditer(text):
        name, body = m.group(1).upper(), m.group(2)
        entries: dict[str, Any] = {}
        # split "KEY = v1, v2, KEY2 = v3" on assignments
        pieces = re.split(r"(\b\w+\s*=)", body)
        key = None
        for piece in pieces:
            if re.fullmatch(r"\s*\w+\s*=", piece or ""):
                key = piece.strip().rstrip("=").strip().upper()
                entries[key] = []
            elif key is not None:
                entries[key].extend(_fortran_scalar(v) for v in _split_values(piece))
        groups[name] = {k: (v[0] if len(v) == 1 else v) for k, v in entries.items()}
    return groups


def case_from_control_in(path, mesh_dir: str | None = None) -> Case:
    """Build a :class:`Case` from a Fortran ``control.in``."""
    g = parse_control_in(path)

    def get(group, key, default=None):
        return g.get(group, {}).get(key, default)

    def aslist(v):
        return list(v) if isinstance(v, list) else [v]

    nbc = int(get("N_BC", "NBC", 1))
    names = aslist(get("BC", "BC_NAME", []))
    phyid = aslist(get("BC", "BC_PHYID", []))
    typ = aslist(get("BC", "BC_TYP", []))
    temp = aslist(get("BC", "BC_TEMP", []))
    xoff = aslist(get("BC", "BC_XOFF", [0.0] * nbc))
    yoff = aslist(get("BC", "BC_YOFF", [0.0] * nbc))

    boundaries = [
        Boundary(name=str(names[i]), phyid=int(phyid[i]),
                 type=BC_NAME_BY_CODE[int(typ[i])],
                 temp=float(temp[i]),
                 xoff=float(xoff[i]) if i < len(xoff) else 0.0,
                 yoff=float(yoff[i]) if i < len(yoff) else 0.0)
        for i in range(nbc)
    ]

    mesh_file = str(get("FILENAME", "FNAME_MSH", "Gmsh.msh"))
    if mesh_dir is not None:
        mesh_file = str((Path(mesh_dir) / Path(mesh_file).name).resolve())

    return Case(
        iteration=Iteration(tol=float(get("ITERATION", "TOL", 1e-8)),
                            tmax=int(get("ITERATION", "TMAX", 10000))),
        scheme=Scheme(accflag=int(get("GSIS", "ACCFLAG", 0))),
        velmesh=VelMesh(npole=int(get("VELMSH", "NPOLE", 20)),
                        nazim=int(get("VELMSH", "NAZIM", 12))),
        dg=DG(deg=int(get("DG", "DEG", 1))),
        flow=Flow(cv=float(get("FLOW", "CV", 0.0)),
                  vg=float(get("FLOW", "VG", 1.0)),
                  tau_r=float(get("FLOW", "TAU_R", 0.0)),
                  tau_n=float(get("FLOW", "TAU_N", 0.0)),
                  tau_thr=float(get("FLOW", "TAU_THR", 1.0))),
        mesh=MeshCfg(file=mesh_file),
        boundaries=boundaries,
        basedir=str(Path(path).parent),
    )
