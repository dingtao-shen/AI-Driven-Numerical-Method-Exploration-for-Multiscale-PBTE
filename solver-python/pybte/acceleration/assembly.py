"""HDG local operators and the global trace matrix -- ``Init_Acceleration_New``.

The macroscopic system carries seven fields per element DOF,

    UQ = [ T, qx, qy, Lxx, Lxy, Lyx, Lyy ]

and three trace fields per face DOF,

    U_TRACE = [ T_hat, qx_hat, qy_hat ]

The element (``local``) problem expresses ``UQ`` in terms of the traces on its
three faces and a source built from the kinetic solution::

    UQ_I = inv(AA_SOL_I) @ ( AA_SRC_I + sum_IL AA_TRACE_{I,IL} @ U_TRACE_face )

Substituting that into the face (``global``) problem gives a sparse system in
the traces alone.  Its matrix depends only on geometry, ``TAU_R``/``TAU_C``
and the stabilisation, so it is assembled and factorised **once**; only the
right-hand side changes from iteration to iteration.

Two deviations from the Fortran are implemented here:

§7.2  Both rescalings are available.  Variant A (the Makefile's) keeps
      ``OO/TAU_R`` in the momentum block, which blows up as ``TAU_R -> 0``;
      variant B multiplies the temperature and stress couplings by ``TAU_R``
      instead, leaving ``OO + ST*PP``.  Algebraically the same system,
      very different conditioning.

§7.4  The Fortran zeroes a dense ``(3*NDOF_FC, N_FCS*3*NDOF_FC)`` row-block
      workspace inside the loop over faces, which makes assembly quadratic in
      face count.  Here each face's row block is accumulated into a small dict
      of column blocks and emitted straight to COO triplets, so assembly is
      linear.
"""
from __future__ import annotations

import numpy as np

from ..constants import BC_NONTHERMALISING, BC_PERIODIC, BC_THERMALISING

__all__ = ["build_local_operators", "build_ba_sol", "assemble_global_matrix"]


def build_local_operators(mesh, integrals, ndof_tri, ndof_fc, cv, tau_r, tau_c,
                          st, variant="A"):
    """Return ``(inv_aa_sol, aa_trace)``.

    ``inv_aa_sol`` : (n_tris, 7*nd, 7*nd)
    ``aa_trace``   : (n_tris, 3, 7*nd, 3*nf)  already premultiplied by the
                     inverse, exactly as the Fortran leaves it
    """
    nd, nf = ndof_tri, ndof_fc
    n_tris = mesh.n_tris
    scale = tau_r if variant == "B" else 1.0

    aa = integrals.int_tri_tri_x                       # (I, M, L)
    bb = integrals.int_tri_tri_y
    cc = np.transpose(aa, (0, 2, 1)).copy()            # CC(M,L) = AA(L,M)
    dd = np.transpose(bb, (0, 2, 1)).copy()
    pp = integrals.int_tri_fc.sum(axis=1)              # sum over the 3 edges
    oo = np.transpose(integrals.int_tri_tri, (2, 0, 1)).copy()   # (I, M, L)

    A = np.zeros((n_tris, 7 * nd, 7 * nd))

    def blk(r, c):
        return A[:, r * nd:(r + 1) * nd, c * nd:(c + 1) * nd]

    # T (energy)
    blk(0, 0)[:] = st[0] * pp
    blk(0, 1)[:] = cc
    blk(0, 2)[:] = dd

    # qx (momentum)
    blk(1, 0)[:] = -aa * cv / 3.0 * scale
    blk(1, 1)[:] = (oo + st[1] * pp) if variant == "B" else (st[1] * pp + oo / tau_r)
    blk(1, 3)[:] = -4.0 / 3.0 * cc
    blk(1, 4)[:] = -dd
    blk(1, 5)[:] = -dd
    blk(1, 6)[:] = 2.0 / 3.0 * cc

    # qy
    blk(2, 0)[:] = -bb * cv / 3.0 * scale
    blk(2, 2)[:] = (oo + st[2] * pp) if variant == "B" else (st[2] * pp + oo / tau_r)
    blk(2, 3)[:] = 2.0 / 3.0 * dd
    blk(2, 4)[:] = -cc
    blk(2, 5)[:] = -cc
    blk(2, 6)[:] = -4.0 / 3.0 * dd

    # stress closure
    tc5 = tau_c * scale / 5.0
    blk(3, 1)[:] = aa * tc5
    blk(3, 3)[:] = oo
    blk(4, 1)[:] = bb * tc5
    blk(4, 4)[:] = oo
    blk(5, 2)[:] = aa * tc5
    blk(5, 5)[:] = oo
    blk(6, 2)[:] = bb * tc5
    blk(6, 6)[:] = oo

    inv_aa_sol = np.linalg.inv(A)

    # ---- AA_TRACE ------------------------------------------------------
    ww = np.transpose(integrals.int_tri_fc_fc, (2, 3, 0, 1)).copy()  # (I, IL, nd, nf)
    nx = mesh.tri_normal[:, :, 0]
    ny = mesh.tri_normal[:, :, 1]

    T = np.zeros((n_tris, 3, 7 * nd, 3 * nf))

    def tblk(r, c):
        return T[:, :, r * nd:(r + 1) * nd, c * nf:(c + 1) * nf]

    nxe = nx[:, :, None, None]
    nye = ny[:, :, None, None]
    tblk(0, 0)[:] = st[0] * ww
    tblk(1, 0)[:] = -nxe * ww * cv / 3.0 * scale
    tblk(1, 1)[:] = st[1] * ww
    tblk(2, 0)[:] = -nye * ww * cv / 3.0 * scale
    tblk(2, 2)[:] = st[2] * ww
    tblk(3, 1)[:] = nxe * ww * tc5
    tblk(4, 1)[:] = nye * ww * tc5
    tblk(5, 2)[:] = nxe * ww * tc5
    tblk(6, 2)[:] = nye * ww * tc5

    aa_trace = np.einsum("ipq,ilqr->ilpr", inv_aa_sol, T, optimize=True)
    return inv_aa_sol, aa_trace


def build_ba_sol(mesh, integrals, ndof_tri, ndof_fc, st, bc_type, variant="A"):
    """Return ``BA_SOL`` with shape ``(n_tris, 3, 3*nf, 7*nd)``.

    This is the face-equation operator: it turns the element fields into the
    numerical trace flux.  Interior and periodic faces get the standard HDG
    average-plus-jump form; a thermalising wall gets the modified rows that
    impose the wall momentum balance, and its temperature rows stay zero
    because the trace temperature is imposed directly from the kinetic
    solution (see ``global_solve``).
    """
    nd, nf = ndof_tri, ndof_fc
    n_tris = mesh.n_tris
    B = np.zeros((n_tris, 3, 3 * nf, 7 * nd))
    bb2 = np.transpose(integrals.int_tri_fc_fc, (2, 3, 1, 0)).copy()  # (I, IL, nf, nd)
    nx = mesh.tri_normal[:, :, 0][:, :, None, None]
    ny = mesh.tri_normal[:, :, 1][:, :, None, None]

    def blk(r, c):
        return B[:, :, r * nf:(r + 1) * nf, c * nd:(c + 1) * nd]

    face_bc = mesh.face_bc[mesh.tri_faces]                     # (n_tris, 3)
    code = np.where(face_bc >= 0, bc_type[np.clip(face_bc, 0, None)], 0)
    interior = (face_bc < 0)
    periodic = (~interior) & (code == BC_PERIODIC)
    wall_th = (~interior) & (code == BC_THERMALISING)
    std = (interior | periodic)[:, :, None, None]

    blk(0, 0)[:] = np.where(std, 0.5 * bb2, 0.0)
    blk(0, 1)[:] = np.where(std, nx / 2.0 / st[0] * bb2, 0.0)
    blk(0, 2)[:] = np.where(std, ny / 2.0 / st[0] * bb2, 0.0)

    blk(1, 1)[:] = np.where(std, 0.5 * bb2, 0.0)
    blk(1, 3)[:] = np.where(std, -2.0 * nx / 3.0 / st[1] * bb2, 0.0)
    blk(1, 4)[:] = np.where(std, -ny / 2.0 / st[1] * bb2, 0.0)
    blk(1, 5)[:] = np.where(std, -ny / 2.0 / st[1] * bb2, 0.0)
    blk(1, 6)[:] = np.where(std, nx / 3.0 / st[1] * bb2, 0.0)

    blk(2, 2)[:] = np.where(std, 0.5 * bb2, 0.0)
    blk(2, 3)[:] = np.where(std, ny / 3.0 / st[2] * bb2, 0.0)
    blk(2, 4)[:] = np.where(std, -nx / 2.0 / st[2] * bb2, 0.0)
    blk(2, 5)[:] = np.where(std, -nx / 2.0 / st[2] * bb2, 0.0)
    blk(2, 6)[:] = np.where(std, -2.0 * ny / 3.0 / st[2] * bb2, 0.0)

    if variant == "A":
        w = wall_th[:, :, None, None]
        blk(1, 1)[:] = np.where(w, bb2, blk(1, 1))
        blk(1, 3)[:] = np.where(w, -4.0 * nx / 3.0 / st[1] * bb2, blk(1, 3))
        blk(1, 4)[:] = np.where(w, -ny / st[1] * bb2, blk(1, 4))
        blk(1, 5)[:] = np.where(w, -ny / st[1] * bb2, blk(1, 5))
        blk(1, 6)[:] = np.where(w, 2.0 * nx / 3.0 / st[1] * bb2, blk(1, 6))

        blk(2, 2)[:] = np.where(w, bb2, blk(2, 2))
        blk(2, 3)[:] = np.where(w, 2.0 * ny / 3.0 / st[2] * bb2, blk(2, 3))
        blk(2, 4)[:] = np.where(w, -nx / st[2] * bb2, blk(2, 4))
        blk(2, 5)[:] = np.where(w, -nx / st[2] * bb2, blk(2, 5))
        blk(2, 6)[:] = np.where(w, -4.0 * ny / 3.0 / st[2] * bb2, blk(2, 6))
    return B


def mirror_permutation(ndof_fc: int) -> np.ndarray:
    """Row permutation that mirrors a face's trace DOFs, per field block.

    Face DOF ``k`` sits a fraction ``k/(NDOF_FC-1)`` of the way from the
    face's node 1 to its node 2.  A periodic partner traversed in the
    opposite sense has its DOF ``k`` at the *other* end, so contributions
    coming across the pair must be reindexed ``k -> NDOF_FC-1-k`` within each
    of the three field blocks.
    """
    base = np.arange(ndof_fc)[::-1]
    return np.concatenate([base + b * ndof_fc for b in range(3)])


def assemble_global_matrix(mesh, integrals, ba_sol, aa_trace, ndof_fc,
                           bc_type, drop_tol=1e-30, periodic_flip=None):
    """Assemble the sparse trace system directly into COO triplets (§7.4).

    Returns ``(csr_matrix, nnz_before_drop)``.
    """
    from scipy.sparse import coo_matrix

    nf3 = 3 * ndof_fc
    n_faces = mesh.n_faces
    n = n_faces * nf3
    perm = mirror_permutation(ndof_fc)

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    vals: list[np.ndarray] = []
    rr, cc = np.meshgrid(np.arange(nf3), np.arange(nf3), indexing="ij")

    aa_tmp = np.zeros((nf3, nf3))
    for f in range(n_faces):
        aa = integrals.int_fc_fc[:, :, f]
        aa_tmp[:] = 0.0
        for k in range(3):
            aa_tmp[k * ndof_fc:(k + 1) * ndof_fc,
                   k * ndof_fc:(k + 1) * ndof_fc] = aa

        blocks: dict[int, np.ndarray] = {f: aa_tmp.copy()}

        bc = int(mesh.face_bc[f])
        flip = bool(periodic_flip[f]) if periodic_flip is not None else False
        if bc < 0:
            t1, t2 = mesh.face_tri[f]
            l1, l2 = mesh.face_lfc[f]
            pairs = [(int(t1), int(l1), False), (int(t2), int(l2), False)]
        elif bc_type[bc] == BC_PERIODIC:
            t1, l1 = mesh.adjacent(f)
            t2, l2 = mesh.adjacent(int(mesh.face_pair[f]))
            # the partner's equation rows are written in *its* face basis, so
            # they must be mirrored back into this face's basis (§7.1)
            pairs = [(t1, l1, False), (t2, l2, flip)]
        else:
            t1, l1 = mesh.adjacent(f)
            pairs = [(t1, l1, False)]

        for tri, lfc, mirror in pairs:
            for il in range(3):
                g = int(mesh.tri_faces[tri, il])
                d = ba_sol[tri, lfc] @ aa_trace[tri, il]
                if mirror:
                    d = d[perm]
                if g in blocks:
                    blocks[g] -= d
                else:
                    blocks[g] = -d

        base_r = f * nf3
        for g, blk in blocks.items():
            keep = np.abs(blk) > drop_tol
            if not keep.any():
                continue
            rows.append(base_r + rr[keep])
            cols.append(g * nf3 + cc[keep])
            vals.append(blk[keep])

    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    vals = np.concatenate(vals)
    K = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    K.sort_indices()
    return K, int(vals.size)
