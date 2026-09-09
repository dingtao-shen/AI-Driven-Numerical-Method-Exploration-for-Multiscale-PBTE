"""Boundary conditions -- ``Boundary_Conditions.f90`` with the dispatch wired up.

In the Fortran, ``Solvers.f90``'s boundary branch is entirely commented out:
the *thermalising* formula is applied unconditionally to every boundary face,
whatever ``BC_TYP`` says.  The diffusely-reflecting machinery exists in
``Calculate_FLUX_WALL`` but is never wired in, and the periodic branch is dead
even though ``Spatial_Mesh.f90`` does the face pairing.

Here all three are implemented and dispatched on ``BC_TYP``:

``1  thermalising``      inflow ``f_w = Cv*T_wall/(4 pi)``.  Identical to the
                         Fortran for every shipped case, since all of their
                         boundaries are type 1.
``2  non-thermalising``  diffusely reflecting and adiabatic: the wall emits
                         the isotropic ``f_w`` that makes the net normal
                         energy flux vanish.  Evaluated from the previous
                         iterate, as ``Calculate_FLUX_WALL`` would have been.
``3  periodic``          inflow from the element behind the paired face; the
                         geometric offset is already baked into
                         ``INT_NODFUNC_TRI_TRI_FC`` at setup.
``4  symmetry``          mesh bookkeeping only in the reference; rejected here
                         rather than silently treated as a hot wall.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._kernels import wall_flux_kernel
from .constants import (BC_NONTHERMALISING, BC_PERIODIC, BC_SYMMETRY,
                        BC_THERMALISING)

__all__ = ["BoundaryData", "build_boundary_data"]


@dataclass
class BoundaryData:
    bc_type: np.ndarray        # (nbc,) int32
    bc_temp: np.ndarray        # (nbc,)
    periodic_tri: np.ndarray   # (n_faces,) int32; element behind the pair, -1
    periodic_flip: np.ndarray  # (n_faces,) bool; pair traversed in reverse
    flux_wall: np.ndarray      # (n_faces, ndof_tri) diffuse wall emission
    has_nonthermalising: bool
    has_periodic: bool

    def update_wall_flux(self, vdf, cxv, cyv, domega, int_tri_fc, mesh) -> None:
        """Refresh the adiabatic-wall emission from the current VDF."""
        if not self.has_nonthermalising:
            return
        wall_flux_kernel(vdf, cxv, cyv, domega,
                         np.ascontiguousarray(int_tri_fc),
                         np.ascontiguousarray(mesh.tri_normal[:, :, 0]),
                         np.ascontiguousarray(mesh.tri_normal[:, :, 1]),
                         mesh.face_bc.astype(np.int32),
                         mesh.face_tri.astype(np.int32),
                         mesh.face_lfc.astype(np.int32),
                         self.bc_type, self.flux_wall)


def boundary_heat_flux(solver) -> np.ndarray:
    """Net outward normal energy flux through every boundary face.

    Uses the same upwind numerical flux the DG scheme itself uses: the
    outgoing half-space is integrated from the element's own trace, the
    incoming half-space from whatever the wall emits.  For a diffusely
    reflecting (adiabatic) wall the emission is constructed precisely so that
    the two cancel, so this is a sharp test of the dispatch, not a
    loose one -- it should return zero to round-off, not to discretisation
    error.

    The diffuse emission is recomputed from the *current* distribution before
    measuring.  The solver caches it from the start of the last step, one
    sweep stale, which on a converged run shows up as a spurious 1e-6
    imbalance -- an artefact of when the cache was filled, not of the boundary
    condition.

    Returns
    -------
    (n_faces,) array; entries for interior faces are zero.
    """
    m = solver.mesh
    ctx = solver.ctx
    cv = solver.case.flow.cv
    from .constants import PI

    fcm = solver.integrals.int_tri_fc                # (I, IL, M, L)
    out = np.zeros(m.n_faces)
    bd = solver.bcdata
    bd.update_wall_flux(solver.vdf, ctx.cxv, ctx.cyv, ctx.domega, fcm, m)

    for f in range(m.n_faces):
        bc = int(m.face_bc[f])
        if bc < 0:
            continue
        tri, lfc = m.adjacent(f)
        nx, ny = m.tri_normal[tri, lfc]
        u = ctx.cxv * nx + ctx.cyv * ny            # c . n_out, per direction
        phi_int = fcm[tri, lfc].sum(axis=0)        # (ndof,) = oint phi_L ds
        trace = solver.vdf[:, tri, :] @ phi_int    # (ndir,) = oint f ds

        outgoing = u > 0.0
        out[f] = float(np.sum(u[outgoing] * ctx.domega[outgoing] * trace[outgoing]))

        typ = bd.bc_type[bc]
        if typ == BC_NONTHERMALISING:
            # FLUX_WALL stores -f_w, not f_w: its numerator is accumulated
            # against the *inward* normal, which flips the sign.  The sweep
            # restores it with ``FW = -FLUX_WALL`` and so must this.
            emitted = -float(bd.flux_wall[f].sum())         # oint f_w ds
        elif typ == BC_THERMALISING:
            emitted = cv / 4.0 / PI * bd.bc_temp[bc] * float(m.face_len[f])
        else:                                                # periodic
            continue
        incoming = ~outgoing
        out[f] += float(np.sum(u[incoming] * ctx.domega[incoming])) * emitted
    return out


def build_boundary_data(mesh, boundaries, ndof_tri: int) -> BoundaryData:
    bc_type = np.array([b.code for b in boundaries], dtype=np.int32)
    bc_temp = np.array([float(b.temp) for b in boundaries], dtype=np.float64)
    bc_xoff = np.array([float(b.xoff) for b in boundaries], dtype=np.float64)
    bc_yoff = np.array([float(b.yoff) for b in boundaries], dtype=np.float64)

    if np.any(bc_type == BC_SYMMETRY):
        raise NotImplementedError(
            "BC_TYP=4 (symmetry) is carried through the mesh reader by the "
            "reference but never implemented in the transport solve; refusing "
            "to treat it as a thermalising wall. Use a periodic pair instead.")

    n_faces = mesh.n_faces
    periodic_tri = np.full(n_faces, -1, dtype=np.int32)
    periodic_flip = np.zeros(n_faces, dtype=bool)
    has_periodic = bool(np.any(bc_type == BC_PERIODIC))
    if has_periodic:
        for f in range(n_faces):
            bc = int(mesh.face_bc[f])
            if bc < 0 or bc_type[bc] != BC_PERIODIC:
                continue
            pair = int(mesh.face_pair[f])
            if pair < 0:
                raise ValueError(f"periodic face {f} has no partner")
            tri, _ = mesh.adjacent(pair)
            periodic_tri[f] = tri

            # A face's 1-D trace basis is parameterised from its own node 1 to
            # its own node 2.  A periodic pair is normally traversed in
            # opposite senses (the two boundaries are walked the same way
            # around the domain), so the two faces' DOF k refer to *mirrored*
            # points.  Record it: the HDG global assembly has to permute the
            # partner's rows, and the reference -- where the periodic branch
            # is dead code -- does not.
            a = mesh.nodes[mesh.face_nodes[f, 0]] + np.array([bc_xoff[bc], bc_yoff[bc]])
            b1 = mesh.nodes[mesh.face_nodes[pair, 0]]
            b2 = mesh.nodes[mesh.face_nodes[pair, 1]]
            if np.allclose(a, b2) and not np.allclose(a, b1):
                periodic_flip[f] = True
            elif not np.allclose(a, b1):
                raise ValueError(
                    f"periodic faces {f} and {pair} do not line up under the "
                    f"declared offset ({bc_xoff[bc]}, {bc_yoff[bc]})")

    flux_wall = np.zeros((n_faces, ndof_tri))
    has_non = bool(np.any(bc_type == BC_NONTHERMALISING))

    return BoundaryData(bc_type=bc_type, bc_temp=bc_temp,
                        periodic_tri=periodic_tri, periodic_flip=periodic_flip,
                        flux_wall=flux_wall,
                        has_nonthermalising=has_non, has_periodic=has_periodic)
