"""The transport sweep -- the port of ``Solvers.f90::DG_Solver_VDF``.

For each direction the elements are visited in ``TRI_ORDER`` and a dense
``NDOF_TRI x NDOF_TRI`` system is solved per element.  The inflow term reads
``VDF`` of the upwind neighbour *after* that neighbour has already been
overwritten in this same pass -- the sweep is Gauss-Seidel in the sweep
ordering, not Jacobi.  Reproducing that in-place semantics is mandatory: a
Jacobi port converges at a different rate and silently invalidates every
iteration count the benchmark grades on.

Two evaluation modes, numerically identical:

``precomputed``  ``A_SOL`` depends only on geometry, direction and ``TAU_C``,
                 all iteration-invariant, so it is LU-factorised once at setup
                 and the sweep is reduced to a triangular solve.  Costs
                 ``N_TRIS * NDIR * NDOF^2 * 8`` bytes (128 MB for the shipped
                 case at DEG=3).
``onthefly``     rebuild and re-factorise every element every iteration, as
                 the Fortran does.  Slower, constant memory.

``sweep_reference.py`` holds a third, pure-numpy/LAPACK implementation used
only to check the kernels.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import _kernels
from ._numba import HAVE_NUMBA, set_num_threads
from .constants import PI

__all__ = ["SweepContext"]


@dataclass
class SweepContext:
    """Contiguous arrays plus the factorised operators, ready for the kernels."""
    order: np.ndarray
    cxv: np.ndarray
    cyv: np.ndarray
    domega: np.ndarray
    cv: float
    vg: float
    tau_r: float
    tau_n: float
    tau_c: float
    mass: np.ndarray
    gx: np.ndarray
    gy: np.ndarray
    fcm: np.ndarray
    ttfc: np.ndarray
    nx: np.ndarray
    ny: np.ndarray
    tri_faces: np.ndarray
    face_bc: np.ndarray
    neighbour: np.ndarray
    periodic_tri: np.ndarray
    bc_type: np.ndarray
    bc_temp: np.ndarray
    flux_wall: np.ndarray
    mode: str = "precomputed"
    kernel: str = "numba"
    lu: np.ndarray | None = field(default=None, repr=False)
    piv: np.ndarray | None = field(default=None, repr=False)
    sweep_count: int = 0

    # -- construction ------------------------------------------------------
    @classmethod
    def build(cls, case, mesh, integrals, vel, sweep_order, bcdata,
              precompute: bool | None = None, max_bytes: float = 8e9):
        cxv, cyv, domega = vel.flat()
        ctx = cls(
            order=np.ascontiguousarray(sweep_order.order, dtype=np.int32),
            cxv=np.ascontiguousarray(cxv), cyv=np.ascontiguousarray(cyv),
            domega=np.ascontiguousarray(domega),
            cv=float(case.flow.cv), vg=float(case.flow.vg),
            tau_r=float(case.flow.tau_r), tau_n=float(case.flow.tau_n),
            tau_c=float(case.tau_c),
            mass=np.ascontiguousarray(integrals.int_tri_tri),
            gx=np.ascontiguousarray(integrals.int_tri_tri_x),
            gy=np.ascontiguousarray(integrals.int_tri_tri_y),
            fcm=np.ascontiguousarray(integrals.int_tri_fc),
            ttfc=np.ascontiguousarray(integrals.int_tri_tri_fc),
            nx=np.ascontiguousarray(mesh.tri_normal[:, :, 0]),
            ny=np.ascontiguousarray(mesh.tri_normal[:, :, 1]),
            tri_faces=np.ascontiguousarray(mesh.tri_faces, dtype=np.int32),
            face_bc=np.ascontiguousarray(mesh.face_bc, dtype=np.int32),
            neighbour=np.ascontiguousarray(mesh.tri_neighbour, dtype=np.int32),
            periodic_tri=bcdata.periodic_tri,
            bc_type=bcdata.bc_type, bc_temp=bcdata.bc_temp,
            flux_wall=bcdata.flux_wall,
        )
        set_num_threads(case.performance.threads)
        ctx.kernel = case.performance.kernel if HAVE_NUMBA else "numpy"

        if precompute is None:
            precompute = case.performance.precompute_inverse
        nbytes = ctx.operator_bytes(mesh.n_tris, case.ndof_tri)
        if precompute and nbytes <= max_bytes and ctx.kernel == "numba":
            ctx.factorise()
        else:
            # the numpy reference kernel builds its own operators inline
            ctx.mode = "onthefly"
        return ctx

    def operator_bytes(self, n_tris: int, ndof: int) -> float:
        return float(self.cxv.size) * n_tris * ndof * ndof * 8.0

    def factorise(self) -> None:
        """LU-factorise ``A_SOL`` for every (direction, element) pair."""
        ndir = self.cxv.size
        nd, _, n_tris = self.mass.shape
        self.lu = np.empty((ndir, n_tris, nd, nd))
        self.piv = np.empty((ndir, n_tris, nd), dtype=np.int32)
        bad = _kernels.build_operators(
            self.order, self.cxv, self.cyv, 1.0 / self.tau_c,
            self.mass, self.gx, self.gy, self.fcm, self.nx, self.ny,
            self.lu, self.piv)
        if bad:
            raise np.linalg.LinAlgError(
                f"{bad} element operators are singular; check TAU_C and the mesh")
        self.mode = "precomputed"

    # -- the sweep ---------------------------------------------------------
    def sweep(self, mom, vdf) -> None:
        """One transport sweep, in place on ``vdf`` (ndir, n_tris, ndof)."""
        if self.kernel == "numpy":
            from .sweep_reference import sweep_reference

            sweep_reference(self, mom, vdf)
            self.sweep_count += 1
            return
        if self.mode == "precomputed":
            _kernels.sweep_precomputed(
                self.order, self.cxv, self.cyv, self.cv, self.vg,
                self.tau_r, self.tau_n, PI,
                self.mass, self.ttfc, self.nx, self.ny, self.tri_faces,
                self.face_bc, self.neighbour, self.periodic_tri,
                self.bc_type, self.bc_temp, self.flux_wall,
                mom.ts, mom.qxs, mom.qys, self.lu, self.piv, vdf)
        else:
            _kernels.sweep_onthefly(
                self.order, self.cxv, self.cyv, self.cv, self.vg,
                self.tau_r, self.tau_n, self.tau_c, PI,
                self.mass, self.gx, self.gy, self.fcm, self.ttfc,
                self.nx, self.ny, self.tri_faces, self.face_bc,
                self.neighbour, self.periodic_tri, self.bc_type,
                self.bc_temp, self.flux_wall, mom.ts, mom.qxs, mom.qys, vdf)
        self.sweep_count += 1

    # -- §7.6 --------------------------------------------------------------
    def true_residual(self, mom, vdf) -> float:
        """``||A_SOL f - A_SRC|| / ||A_SRC||`` in the solid-angle-weighted
        L2 norm, evaluated at the current state with *no* sweep applied.

        This asks one question: do ``vdf`` and the moments in ``mom`` jointly
        satisfy the discrete transport system?  It is zero exactly at the
        fixed point of whichever scheme produced them.

        Read the two schemes differently.

        **CIS.**  The sweep leaves ``A f = b(moments_before)``, and ``mom`` is
        then recomputed from that ``f``, so what this measures is
        ``b(moments_before) - b(moments_after)`` -- the change in the moments,
        which is what ``residual_iterate`` already measures.  The two agree to
        within a normalisation (1.52e-06 vs 1.52e-06 on the truncated
        ``tau_R = 1e-3`` run).  **Neither reveals the error**: see
        ``RunRecord.error_estimate`` for that.

        **GSIS.**  The correction moves ``vdf`` and the moments together to a
        state that does *not* satisfy the transport system, and this is what
        detects it -- 3.2e-03 at ``tau_R = 1e-1`` where CIS reaches 4.3e-14.
        That is the evidence behind docs/FORTRAN_ISSUES.md #5.
        """
        ndir = self.cxv.size
        num = np.zeros(ndir)
        den = np.zeros(ndir)
        _kernels.transport_residual(
            self.cxv, self.cyv, self.domega, self.cv, self.vg,
            self.tau_r, self.tau_n, self.tau_c, PI,
            self.mass, self.gx, self.gy, self.fcm, self.ttfc,
            self.nx, self.ny, self.tri_faces, self.face_bc, self.neighbour,
            self.periodic_tri, self.bc_type, self.bc_temp, self.flux_wall,
            mom.ts, mom.qxs, mom.qys, vdf, num, den)
        d = float(den.sum())
        if d <= 0.0:
            return float("nan")
        return float(np.sqrt(num.sum() / d))
