"""Global trace problem -- ``Global_Problem_Solver_ACC`` + ``Init_PARDISO``.

The matrix is iteration-invariant, so it is factorised exactly once with
``scipy.sparse.linalg.splu`` and every iteration only calls ``.solve()`` on a
new right-hand side.  That mirrors PARDISO's phase separation (12 then 33)
and is what keeps a GSIS iteration in the same cost bracket as a CIS one.

The right-hand side has two parts:

* the interior/periodic contribution ``BA_SOL @ AA_SRC`` from the one or two
  elements that touch the face;
* on a physical wall, the trace fields are imposed **from the kinetic
  solution**, not from the wall temperature: the face rows are overwritten
  with the DVM's own wall traces.  A thermalising wall imposes ``T_hat``
  only, an adiabatic one imposes all three.  (The Fortran keeps the
  wall-temperature alternative commented out just above; taking the trace
  from the DVM is what makes GSIS preserve the kinetic fixed point instead
  of converging to a Fourier solution.)
"""
from __future__ import annotations

import numpy as np

from .._numba import njit, prange
from ..constants import BC_NONTHERMALISING, BC_PERIODIC, BC_THERMALISING

__all__ = ["GlobalSolver", "wall_trace_moments"]


@njit(cache=True, parallel=True)
def _wall_trace_kernel(vdf, tri_fc_fc, cxv, cyv, domega, face_list,
                       tri_list, lfc_list, nd, out):
    """``out[k, 0:nf] = T_wall``, ``[nf:2nf] = qx_wall``, ``[2nf:3nf] = qy``."""
    nf = out.shape[1] // 3
    ndir = vdf.shape[0]
    for k in prange(face_list.shape[0]):
        tri = tri_list[k]
        lfc = lfc_list[k]
        for m in range(nf):
            tbc = 0.0
            qxbc = 0.0
            qybc = 0.0
            for d in range(ndir):
                t = 0.0
                for l in range(nd):
                    t += vdf[d, tri, l] * tri_fc_fc[l, m, tri, lfc]
                tbc += t * domega[d]
                qxbc += cxv[d] * t * domega[d]
                qybc += cyv[d] * t * domega[d]
            out[k, m] = tbc
            out[k, nf + m] = qxbc
            out[k, 2 * nf + m] = qybc


def wall_trace_moments(vdf, tri_fc_fc, cxv, cyv, domega, face_list, tri_list,
                       lfc_list, ndof_tri, ndof_fc):
    out = np.zeros((face_list.size, 3 * ndof_fc))
    if face_list.size:
        _wall_trace_kernel(vdf, tri_fc_fc, cxv, cyv, domega, face_list,
                           tri_list, lfc_list, ndof_tri, out)
    return out


class GlobalSolver:
    """Owns the assembled matrix and its LU factorisation."""

    def __init__(self, mesh, integrals, ba_sol, ndof_tri, ndof_fc, bc_type,
                 matrix, cv, variant="A", periodic_flip=None):
        from .assembly import mirror_permutation

        self.variant = variant
        self.perm = mirror_permutation(ndof_fc)
        self.periodic_flip = (periodic_flip if periodic_flip is not None
                              else np.zeros(mesh.n_faces, dtype=bool))
        self.mesh = mesh
        self.integrals = integrals
        self.ba_sol = ba_sol
        self.nd = ndof_tri
        self.nf = ndof_fc
        self.nf3 = 3 * ndof_fc
        self.bc_type = bc_type
        self.cv = cv
        self.K = matrix
        self._lu = None
        self.solve_count = 0

        n_faces = mesh.n_faces
        # face -> contributing (element, local face) pairs, in Fortran order
        p1_tri = np.full(n_faces, -1, dtype=np.int64)
        p1_lfc = np.full(n_faces, -1, dtype=np.int64)
        p2_tri = np.full(n_faces, -1, dtype=np.int64)
        p2_lfc = np.full(n_faces, -1, dtype=np.int64)
        for f in range(n_faces):
            bc = int(mesh.face_bc[f])
            if bc < 0:
                p1_tri[f], p2_tri[f] = mesh.face_tri[f]
                p1_lfc[f], p2_lfc[f] = mesh.face_lfc[f]
            elif bc_type[bc] == BC_PERIODIC:
                p1_tri[f], p1_lfc[f] = mesh.adjacent(f)
                p2_tri[f], p2_lfc[f] = mesh.adjacent(int(mesh.face_pair[f]))
            else:
                p1_tri[f], p1_lfc[f] = mesh.adjacent(f)
        self.p1_tri, self.p1_lfc = p1_tri, p1_lfc
        self.p2_tri, self.p2_lfc = p2_tri, p2_lfc
        self.has_p2 = p2_tri >= 0

        bc_of_face = mesh.face_bc
        codes = np.where(bc_of_face >= 0, bc_type[np.clip(bc_of_face, 0, None)], 0)
        self.wall_th = np.flatnonzero((bc_of_face >= 0) & (codes == BC_THERMALISING))
        self.wall_ad = np.flatnonzero((bc_of_face >= 0) & (codes == BC_NONTHERMALISING))
        self.wall_all = np.concatenate([self.wall_th, self.wall_ad])
        self.wall_tri = p1_tri[self.wall_all]
        self.wall_lfc = p1_lfc[self.wall_all]
        self.n_th = self.wall_th.size
        self.wall_nx = mesh.tri_normal[self.wall_tri, self.wall_lfc, 0]
        self.wall_ny = mesh.tri_normal[self.wall_tri, self.wall_lfc, 1]

    # -- factorisation -----------------------------------------------------
    def factorise(self):
        from scipy.sparse.linalg import splu

        self._lu = splu(self.K.tocsc())
        return self

    @property
    def nnz(self) -> int:
        return int(self.K.nnz)

    # -- right-hand side ---------------------------------------------------
    def rhs(self, aa_src, vdf, cxv, cyv, domega):
        nf3 = self.nf3
        n_faces = self.mesh.n_faces
        ffa = np.zeros((nf3, n_faces))

        c1 = np.einsum("fpq,qf->pf",
                       self.ba_sol[self.p1_tri, self.p1_lfc],
                       aa_src[:, self.p1_tri], optimize=True)
        ffa[:] = c1
        sel = self.has_p2
        if sel.any():
            c2 = np.einsum("fpq,qf->pf",
                           self.ba_sol[self.p2_tri[sel], self.p2_lfc[sel]],
                           aa_src[:, self.p2_tri[sel]], optimize=True)
            # mirror the rows that arrive across a reversed periodic pair
            fl = self.periodic_flip[sel]
            if fl.any():
                c2[:, fl] = c2[self.perm][:, fl]
            ffa[:, sel] += c2

        if self.wall_all.size:
            mom = wall_trace_moments(vdf, self.integrals.int_tri_fc_fc,
                                     cxv, cyv, domega, self.wall_all,
                                     self.wall_tri, self.wall_lfc,
                                     self.nd, self.nf)
            nf = self.nf
            k_th = slice(0, self.n_th)
            k_ad = slice(self.n_th, None)
            th = self.wall_all[k_th]
            ad = self.wall_all[k_ad]

            # Thermalising wall.  Variant A imposes only the trace
            # temperature -- its BA_SOL keeps the wall momentum-balance rows.
            # Variant B has no wall rows in BA_SOL at all, so all three trace
            # fields have to be imposed from the kinetic solution.
            if th.size:
                ffa[0:nf, th] = mom[k_th, 0:nf].T / self.cv
                if self.variant == "B":
                    ffa[nf:2 * nf, th] = mom[k_th, nf:2 * nf].T
                    ffa[2 * nf:3 * nf, th] = mom[k_th, 2 * nf:3 * nf].T

            # Adiabatic wall.  Variant B projects the imposed heat flux onto
            # the wall tangent (n.q_hat = 0 exactly); variant A imposes the
            # unrotated DVM flux.
            if ad.size:
                ffa[0:nf, ad] = mom[k_ad, 0:nf].T / self.cv
                qxb = mom[k_ad, nf:2 * nf].T
                qyb = mom[k_ad, 2 * nf:3 * nf].T
                if self.variant == "B":
                    nx = self.wall_nx[k_ad]
                    ny = self.wall_ny[k_ad]
                    ffa[nf:2 * nf, ad] = qxb * (ny * ny) - qyb * (nx * ny)
                    ffa[2 * nf:3 * nf, ad] = -qxb * (nx * ny) + qyb * (nx * nx)
                else:
                    ffa[nf:2 * nf, ad] = qxb
                    ffa[2 * nf:3 * nf, ad] = qyb
        return ffa

    # -- solve -------------------------------------------------------------
    def solve(self, ffa):
        if self._lu is None:
            self.factorise()
        x = self._lu.solve(ffa.ravel(order="F"))
        self.solve_count += 1
        return x.reshape(self.nf3, self.mesh.n_faces, order="F")
