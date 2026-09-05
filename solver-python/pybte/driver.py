"""The outer iteration -- the port of ``Callaway_2D_2V_DG.f90``.

::

    STEP = 1
    DO
        DG_Solver_VDF                       ! transport sweep
        IF CIS  : Calculate_Macro_Properties
        IF GSIS : HoT source; global trace solve; local recovery; correction
        Calculate_Residual_T
        EXIT when RESIDUAL < TOL or STEP >= TMAX
    END DO

Nothing about that structure changes in the port.  What is added around it:
an explicit restart contract (§7.3), an optional true residual (§7.6), and a
:class:`~pybte.io_output.RunRecord` complete enough that a verifier never has
to parse stdout (§8).
"""
from __future__ import annotations

import hashlib
import time
import warnings
from pathlib import Path

import numpy as np

from .basis import build_basis
from .bc import build_boundary_data
from .config import Case
from .integration import build_integrals
from .io_output import (RunRecord, append_runtime_log, locate_points,
                        sample_fields, sampling_grid, write_tecplot)
from .mesh import build_mesh
from .moments import Moments, compute_moments, residual_iterate
from .ordering import build_sweep_order, validate_sweep_order
from .sweep import SweepContext
from .velocity import build_velocity_mesh

__all__ = ["Solver"]


class Solver:
    """Set up a case and run it to convergence.

    >>> case = Case.from_yaml("cases/cavity_tauR1e-3_cis.yaml")
    >>> record = Solver(case).run()
    >>> record.iterations, record.converged
    """

    def __init__(self, case: Case, validate_order: bool = False,
                 on_cycle: str | None = None):
        t0 = time.perf_counter()
        self.case = case
        c = case
        on_cycle = on_cycle or c.scheme.on_cycle

        if not c.boundaries:
            raise ValueError("the case declares no boundaries")

        self.mesh = build_mesh(c.mesh_path(),
                               [b.phyid for b in c.boundaries],
                               c.bc_codes(),
                               [b.name for b in c.boundaries],
                               [b.xoff for b in c.boundaries],
                               [b.yoff for b in c.boundaries])
        self.basis = build_basis(c.dg.deg)
        self.vel = build_velocity_mesh(c.velmesh.npole, c.velmesh.nazim, c.flow.vg)
        self.integrals = build_integrals(self.mesh, self.basis, c.np_tri, c.np_fc,
                                         bc_codes=c.bc_codes(),
                                         bc_xoff=[b.xoff for b in c.boundaries],
                                         bc_yoff=[b.yoff for b in c.boundaries])
        self.order = build_sweep_order(self.mesh, self.vel, on_cycle=on_cycle)
        if validate_order:
            validate_sweep_order(self.mesh, self.vel, self.order)
        self.bcdata = build_boundary_data(self.mesh, c.boundaries, c.ndof_tri)

        self.ndir = self.vel.ndir
        self.n_tris = self.mesh.n_tris
        self.ndof = c.ndof_tri

        self.ctx = SweepContext.build(c, self.mesh, self.integrals, self.vel,
                                      self.order, self.bcdata)
        self.factorisation_count = (self.ndir * self.n_tris
                                    if self.ctx.mode == "precomputed" else 0)

        self.vdf = np.zeros((self.ndir, self.n_tris, self.ndof))
        self.mom = Moments.zeros(self.ndof, self.n_tris)
        self.t_old = np.zeros(self.n_tris)
        self.restart_info = self._apply_restart()

        self.acc = None
        if c.scheme.accflag == 1:
            from .acceleration import Acceleration
            self.acc = Acceleration(self)

        self.setup_seconds = time.perf_counter() - t0

    # -- restart (§7.3) ----------------------------------------------------
    def _restart_stem(self) -> str:
        c = self.case
        return (f"VDFP{c.dg.deg}T{self.n_tris:3d}"
                f"NP{c.velmesh.npole:2d}NA{c.velmesh.nazim:2d}").replace(" ", "")

    def _apply_restart(self) -> dict:
        r = self.case.restart
        default = Path(self.case.basedir) / f"{self._restart_stem()}.out"
        if not r.enabled:
            if r.error_on_stale and default.exists():
                raise FileExistsError(
                    f"a restart file {default} is present but restart is disabled. "
                    "The Fortran would have read it silently and changed the "
                    "iteration count with no warning. Delete it, or set "
                    "restart.enabled or restart.error_on_stale explicitly.")
            return {"enabled": False, "path": None, "sha256": None}

        path = Path(r.path) if r.path else default
        if not path.is_absolute():
            path = Path(self.case.basedir) / path
        if not path.exists():
            raise FileNotFoundError(f"restart enabled but {path} does not exist")
        raw = path.read_bytes()
        want = self.ndof * self.n_tris * self.ndir * 8
        if len(raw) != want:
            raise ValueError(f"{path}: {len(raw)} bytes, expected {want} "
                             "(wrong DEG / mesh / angular resolution?)")
        arr = np.frombuffer(raw, dtype="<f8").reshape(
            (self.ndof, self.n_tris, self.case.velmesh.npole,
             self.case.velmesh.nazim), order="F")
        self.vdf = np.ascontiguousarray(
            np.transpose(arr, (3, 2, 1, 0)).reshape(self.ndir, self.n_tris, self.ndof))
        return {"enabled": True, "path": str(path),
                "sha256": hashlib.sha256(raw).hexdigest()}

    def save_restart(self, path=None) -> Path:
        path = Path(path) if path else Path(self.case.basedir) / f"{self._restart_stem()}.out"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(self.vdf_fortran().ravel(order="F").tobytes())
        return path

    # -- views -------------------------------------------------------------
    def vdf_fortran(self) -> np.ndarray:
        """``VDF(NDOF_TRI, N_TRIS, NPOLE, NAZIM)``, for stage comparison."""
        npole, nazim = self.case.velmesh.npole, self.case.velmesh.nazim
        return np.transpose(
            self.vdf.reshape(nazim, npole, self.n_tris, self.ndof), (3, 2, 1, 0))

    # -- one outer iteration ------------------------------------------------
    def step(self) -> None:
        if self.bcdata.has_nonthermalising:
            self.bcdata.update_wall_flux(self.vdf, self.ctx.cxv, self.ctx.cyv,
                                         self.ctx.domega,
                                         self.integrals.int_tri_fc, self.mesh)
        self.ctx.sweep(self.mom, self.vdf)
        if self.ctx.mode == "onthefly":
            self.factorisation_count += self.ndir * self.n_tris
        if self.acc is None:
            compute_moments(self.vdf, self.ctx.cxv, self.ctx.cyv, self.ctx.domega,
                            self.integrals.int_tri, self.case.flow.cv, self.mom,
                            zero_qy=not self.case.compat.qy_accumulation)
        else:
            self.acc.apply()

    # -- the outer loop -----------------------------------------------------
    def run(self, callback=None, progress_every: int = 0) -> RunRecord:
        c = self.case
        tol, tmax = c.iteration.tol, c.iteration.tmax
        want_true = c.iteration.true_residual
        every = max(1, int(c.iteration.true_residual_every))

        residuals: list[float] = []
        masses: list[float] = []
        true_res: list[float] = []

        t0 = time.perf_counter()
        step = 0
        converged = False
        while True:
            step += 1
            self.step()
            res = residual_iterate(self.mom.temp, self.t_old)
            self.t_old[:] = self.mom.temp
            residuals.append(res)
            masses.append(self.mom.mass)
            if want_true:
                true_res.append(self.ctx.true_residual(self.mom, self.vdf)
                                if (step % every == 0) else np.nan)
            if progress_every and step % progress_every == 0:
                print(f"  step {step:8d}  residual {res:.6e}  mass {self.mom.mass:.10f}",
                      flush=True)
            if callback is not None:
                callback(step, res, self)
            if res < tol:
                converged = True
                break
            if step >= tmax:
                break
        wall = time.perf_counter() - t0

        rec = RunRecord(
            iterations=step, converged=converged,
            residual_history=np.asarray(residuals),
            residual_true=np.asarray(true_res) if want_true else None,
            mass_history=np.asarray(masses),
            temp=self.mom.temp.copy(), qx=self.mom.qx.copy(), qy=self.mom.qy.copy(),
            temp_dofs=self.mom.ts.copy(), qx_dofs=self.mom.qxs.copy(),
            qy_dofs=self.mom.qys.copy(),
            sweep_count=self.ctx.sweep_count,
            factorisation_count=self.factorisation_count,
            wall_clock=wall, setup_wall_clock=self.setup_seconds,
            config_hash=c.config_hash,
            scheme="GSIS" if c.scheme.accflag == 1 else "CIS",
            acc_variant=c.scheme.acc_variant if c.scheme.accflag == 1 else None,
            tol=tol, tmax=tmax, tau_r=c.flow.tau_r, tau_n=c.flow.tau_n,
            tau_c=c.tau_c, deg=c.dg.deg, npole=c.velmesh.npole,
            nazim=c.velmesh.nazim, n_tris=self.n_tris, n_faces=self.mesh.n_faces,
            mesh_file=str(c.mesh_path()),
            restart=self.restart_info,
            environment=RunRecord.environment_fingerprint(),
            diagnostics={
                "sweep_mode": self.ctx.mode,
                "cycles_broken": int(self.order.cycles_broken.sum()),
                "hmin": float(self.mesh.hmin),
                "n_faces_b": int(self.mesh.n_faces_b),
                "boundary_types": [b.type for b in c.boundaries],
            },
        )
        if self.acc is not None:
            rec.diagnostics.update(self.acc.diagnostics())
        return rec

    # -- output -------------------------------------------------------------
    def write_outputs(self, rec: RunRecord) -> dict:
        c = self.case
        out = c.output_path()
        out.mkdir(parents=True, exist_ok=True)
        written = {}

        if c.output.run_record:
            written["run_record"] = str(rec.to_json(out / "run_record.json"))
        if c.output.runtime_log:
            written["runtime_log"] = str(append_runtime_log(
                out / "RunTime.txt", rec.scheme, c.flow.tau_r, c.flow.tau_n,
                self.n_tris, c.velmesh.npole, c.velmesh.nazim,
                rec.iterations, rec.wall_clock))
        if c.output.field:
            px, py = sampling_grid(c.output.field_nx, c.output.field_ny)
            triid = locate_points(self.mesh, px, py)
            missing = int(np.count_nonzero(triid < 0))
            if missing:
                warnings.warn(f"{missing} sampling points fell outside the mesh")
            uq = self.acc.uq if self.acc is not None else None
            pT, pqx, pqy, pxx, pxy, pyy = sample_fields(
                self.mesh, self.basis, self.vdf, self.ctx.cxv, self.ctx.cyv,
                self.ctx.domega, c.flow.cv, px, py, triid, uq=uq,
                ndof_tri=self.ndof)
            if c.output.npz:
                np.savez_compressed(out / "field.npz", px=px, py=py, triid=triid,
                                    T=pT, qx=pqx, qy=pqy,
                                    Nxx=pxx, Nxy=pxy, Nyy=pyy)
                written["field_npz"] = str(out / "field.npz")
            if c.output.tecplot:
                written["field_dat"] = str(write_tecplot(
                    out / "field.dat", px, py, [pT, pqx, pqy, pxx, pxy, pyy],
                    ["T", "qx", "qy", "Nxx", "Nxy", "Nyy"]))
        return written
