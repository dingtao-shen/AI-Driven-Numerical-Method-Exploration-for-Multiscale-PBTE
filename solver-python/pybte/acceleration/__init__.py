"""GSIS -- general synthetic iterative scheme (HDG macroscopic acceleration).

One outer iteration adds four stages after the transport sweep:

1. :func:`~pybte.acceleration.hot_source.hot_source` -- build the
   higher-order-term source from the kinetic solution;
2. :class:`~pybte.acceleration.global_solve.GlobalSolver` -- solve the sparse
   trace system (matrix factorised once at setup);
3. :func:`~pybte.acceleration.local_solve.local_solve` -- recover the seven
   macroscopic fields per element;
4. :func:`~pybte.acceleration.correction.correct_vdf` -- push the macroscopic
   solution back into the distribution, damped by ``beta``.

The whole thing is behind this one class so that the CIS path and the GSIS
path stay cleanly separable at the module boundary (goal G4): with
``accflag = 0`` nothing in this package is imported at all.
"""
from __future__ import annotations

import numpy as np

from ..constants import PI
from .assembly import assemble_global_matrix, build_ba_sol, build_local_operators
from .correction import correct_vdf
from .global_solve import GlobalSolver
from .hot_source import hot_source
from .local_solve import local_solve

__all__ = ["Acceleration", "build_local_operators", "build_ba_sol",
           "assemble_global_matrix", "hot_source", "local_solve",
           "correct_vdf", "GlobalSolver"]


class Acceleration:
    """Everything GSIS needs, built once from a :class:`~pybte.driver.Solver`."""

    def __init__(self, solver):
        c = solver.case
        self.solver = solver
        self.variant = c.scheme.acc_variant
        self.nd = c.ndof_tri
        self.nf = c.ndof_fc

        st = np.asarray(c.scheme.stabilisation, dtype=float).copy()
        if c.scheme.scale_by_hmin:
            st[1] /= solver.mesh.hmin
            st[2] /= solver.mesh.hmin
        self.st = st

        self.inv_aa_sol, self.aa_trace = build_local_operators(
            solver.mesh, solver.integrals, self.nd, self.nf,
            c.flow.cv, c.flow.tau_r, c.tau_c, st, self.variant)

        self.ba_sol = build_ba_sol(solver.mesh, solver.integrals, self.nd,
                                   self.nf, st, solver.bcdata.bc_type,
                                   self.variant)

        drop = 1e-30 if self.variant == "A" else 0.0
        self.K, self.nnz = assemble_global_matrix(
            solver.mesh, solver.integrals, self.ba_sol, self.aa_trace,
            self.nf, solver.bcdata.bc_type, drop_tol=drop,
            periodic_flip=solver.bcdata.periodic_flip)

        self.gsolver = GlobalSolver(solver.mesh, solver.integrals, self.ba_sol,
                                    self.nd, self.nf, solver.bcdata.bc_type,
                                    self.K, c.flow.cv, variant=self.variant,
                                    periodic_flip=solver.bcdata.periodic_flip
                                    ).factorise()

        self.aa_src = np.zeros((7 * self.nd, solver.n_tris))
        self.ffa = np.zeros((3 * self.nf, solver.mesh.n_faces))
        self.u_trace = np.zeros((3 * self.nf, solver.mesh.n_faces))
        self.uq = np.zeros((7 * self.nd, solver.n_tris))
        self._hot_work = None

        self.gx_t = np.ascontiguousarray(solver.integrals.int_tri_tri_x)
        self.gy_t = np.ascontiguousarray(solver.integrals.int_tri_tri_y)
        self.scale = c.flow.tau_r if self.variant == "B" else 1.0
        self.use_hmin = (self.variant == "A")

    # -- one acceleration cycle -------------------------------------------
    def apply(self) -> None:
        s = self.solver
        c = s.case
        self.aa_src, self._hot_work = hot_source(
            s.vdf, self.gx_t, self.gy_t, s.ctx.cxv, s.ctx.cyv, s.ctx.domega,
            c.tau_c, self.nd, self.inv_aa_sol, self.scale, self._hot_work)

        self.ffa = self.gsolver.rhs(self.aa_src, s.vdf, s.ctx.cxv, s.ctx.cyv,
                                    s.ctx.domega)
        self.u_trace = self.gsolver.solve(self.ffa)
        self.uq = local_solve(self.aa_src, self.aa_trace, self.u_trace,
                              s.mesh.tri_faces)
        correct_vdf(s.vdf, self.uq, s.ctx.cxv, s.ctx.cyv, s.ctx.domega,
                    s.integrals.int_tri, s.mesh.tri_hmin,
                    c.flow.cv, c.flow.vg, c.flow.tau_r, c.flow.tau_n,
                    c.tau_c, c.flow.tau_thr, PI, self.nd, s.mom,
                    use_hmin=self.use_hmin)

    # -- diagnostics / cross-validation ------------------------------------
    def diagnostics(self) -> dict:
        return {
            "acc_variant": self.variant,
            "global_matrix_n": int(self.K.shape[0]),
            "global_matrix_nnz": int(self.K.nnz),
            "global_solve_count": int(self.gsolver.solve_count),
            "stabilisation": [float(v) for v in self.st],
        }

    def global_matrix(self):
        return self.K

    def inv_aa_sol_fortran(self) -> np.ndarray:
        return np.transpose(self.inv_aa_sol, (1, 2, 0))

    def aa_trace_fortran(self) -> np.ndarray:
        return np.transpose(self.aa_trace, (2, 3, 1, 0))

    def ba_sol_fortran(self) -> np.ndarray:
        return np.transpose(self.ba_sol, (2, 3, 1, 0))

    def condition_estimate(self) -> float:
        """1-norm condition estimate of the global trace matrix.

        Used to compare the two acceleration variants (§7.2): they are
        algebraically the same system and differ only in conditioning as
        ``TAU_R -> 0``.  Both factors are Higham-Tisseur estimates, so this is
        an estimate of a bound, not the condition number itself.
        """
        from scipy.sparse.linalg import LinearOperator, onenormest

        if self.gsolver._lu is None:
            self.gsolver.factorise()
        lu = self.gsolver._lu
        n = self.K.shape[0]
        inv = LinearOperator((n, n), dtype=float, matvec=lu.solve,
                             rmatvec=lambda b: lu.solve(b, "T"))
        return float(onenormest(self.K) * onenormest(inv))
