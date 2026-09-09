"""Krylov solution of the outer iteration.

One outer step of source iteration is the affine map

    u  ->  F(u) = M( A^-1 b(u) ) = T u + g

on the state ``u``: build the equilibrium source from the moments, invert
the transport operator exactly by one sweep, take moments.  Iterating it is
Richardson's method on

    (I - T) u = g,

which converges at the spectral radius of ``T``.  The Callaway collision
operator conserves energy exactly, so there is no absorption, the scattering
ratio is one, and ``rho(T) = 1 - O((l_mfp / L)^2)``: at ``tau_R = 1e-4`` the
iteration reports a small *step* while sitting nowhere near the answer.

Solving the same system with GMRES costs one sweep per matrix-vector product
and converges on the whole Krylov space instead of the last iterate.  The
fixed point is untouched -- ``(I - T) u = g`` is the equation source
iteration was crawling towards, so its solution *is* the unaccelerated
solver's answer, to round-off rather than to a tolerance.

What the state has to contain
-----------------------------
Everything the sweep reads back from the previous iterate.  With
thermalising walls that is the three moment fields alone: ``TRI_ORDER`` is a
genuine topological order, so the sweep is an exact block-triangular solve
and the post-sweep distribution depends on the moments and the boundary
inflow only.  A diffusely reflecting wall re-emits what it received, so its
emission is appended to the state; it is linear in the distribution and the
map stays affine.  Two things break exactness and are refused by
:func:`applicable`: a direction whose sweep graph had a cycle broken (part
of the inflow is lagged) and periodic faces (the partner element cannot be
placed upwind, so its outflow is lagged too).
"""
from __future__ import annotations

import time

import numpy as np

from .moments import compute_moments

__all__ = ["FixedPointOperator", "applicable", "gmres_full", "run_krylov",
           "MomentPreconditioner"]

_ENOUGH_FACTOR = 1.0e-2     # past this the true residual is reliably under tol
_FLOOR = 1.0e-14            # round-off, relative to ||g||


def applicable(solver) -> tuple[bool, str]:
    """Whether one sweep is an exact evaluation of ``F``.  ``(ok, reason)``."""
    if int(solver.order.cycles_broken.sum()):
        return False, ("the sweep ordering has broken cycles: part of the "
                       "inflow is lagged, so one sweep is not a function of "
                       "the state alone")
    if solver.bcdata.has_periodic:
        return False, ("periodic faces take their inflow from a partner the "
                       "ordering cannot place upwind, so one sweep is not a "
                       "function of the state alone")
    return True, ""


class FixedPointOperator:
    """``u -> F(u)`` and ``v -> (I - T) v``, one sweep apiece."""

    def __init__(self, solver):
        self.s = solver
        self.nd = solver.ndof
        self.nt = solver.n_tris
        self.nmom = 3 * self.nd * self.nt
        self.wall = bool(solver.bcdata.has_nonthermalising)
        self.nwall = int(solver.bcdata.flux_wall.size) if self.wall else 0
        self.n = self.nmom + self.nwall
        self.sweeps = 0

    # -- packing ----------------------------------------------------------
    def scatter(self, u) -> None:
        m = self.s.mom
        k = self.nd * self.nt
        np.copyto(m.ts, u[:k].reshape(m.ts.shape))
        np.copyto(m.qxs, u[k:2 * k].reshape(m.qxs.shape))
        np.copyto(m.qys, u[2 * k:3 * k].reshape(m.qys.shape))
        if self.wall:
            fw = self.s.bcdata.flux_wall
            np.copyto(fw, u[self.nmom:].reshape(fw.shape))

    def gather(self) -> np.ndarray:
        m = self.s.mom
        u = np.empty(self.n)
        k = self.nd * self.nt
        u[:k] = m.ts.reshape(-1)
        u[k:2 * k] = m.qxs.reshape(-1)
        u[2 * k:3 * k] = m.qys.reshape(-1)
        if self.wall:
            u[self.nmom:] = self.s.bcdata.flux_wall.reshape(-1)
        return u

    # -- the map ----------------------------------------------------------
    def apply(self, u) -> np.ndarray:
        """``F(u)``: sweep from the state ``u``, take moments (and emission)."""
        s = self.s
        self.scatter(u)
        s.ctx.sweep(s.mom, s.vdf)
        self.sweeps += 1
        if s.ctx.mode == "onthefly":
            s.factorisation_count += s.ndir * s.n_tris
        compute_moments(s.vdf, s.ctx.cxv, s.ctx.cyv, s.ctx.domega,
                        s.integrals.int_tri, s.case.flow.cv, s.mom)
        if self.wall:
            s.bcdata.update_wall_flux(s.vdf, s.ctx.cxv, s.ctx.cyv, s.ctx.domega,
                                      s.integrals.int_tri_fc, s.mesh)
        return self.gather()

    def state(self) -> np.ndarray:
        """The current state, without sweeping (for a restart)."""
        s = self.s
        compute_moments(s.vdf, s.ctx.cxv, s.ctx.cyv, s.ctx.domega,
                        s.integrals.int_tri, s.case.flow.cv, s.mom)
        if self.wall:
            s.bcdata.update_wall_flux(s.vdf, s.ctx.cxv, s.ctx.cyv, s.ctx.domega,
                                      s.integrals.int_tri_fc, s.mesh)
        return self.gather()


# ---------------------------------------------------------------------------
# GMRES
# ---------------------------------------------------------------------------
def _back_substitute(r, y):
    out = np.array(y, dtype=float)
    for i in range(out.size - 1, -1, -1):
        s = out[i] - r[i, i + 1:] @ out[i + 1:]
        out[i] = s / r[i, i] if r[i, i] != 0.0 else 0.0
    return out


def gmres_full(matvec, b, x0=None, tol_abs=0.0, maxiter=None, max_basis=None,
               precond=None, callback=None):
    """GMRES for ``A x = b`` with ``A`` given as ``matvec``.

    ``max_basis`` bounds the Krylov dimension; when it is at least ``maxiter``
    the method never restarts.  Restarting throws the Krylov space away, and
    on this problem that costs an order of magnitude at small Knudsen number
    (6 343 sweeps with restart 600 against 1 324 without, ``tau_R = 1e-4``),
    so the default is *no* restart with the basis grown by doubling and a
    memory cap deciding where restarts become unavoidable.

    Orthogonalisation is classical Gram-Schmidt applied twice: one pass loses
    orthogonality at the Krylov dimensions this problem reaches (over a
    thousand), two passes match modified Gram-Schmidt and stay BLAS-2.  The
    Hessenberg matrix is kept triangular by Givens rotations so the residual
    of every product is known without forming the iterate.

    ``callback(res)`` is called once per product with the absolute residual;
    returning true stops the solve.  Returns ``(x, converged, n_matvec)``.
    """
    n = b.shape[0]
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float)
    maxiter = int(n if maxiter is None else maxiter)
    m_cap = maxiter if max_basis is None else max(1, min(int(max_basis), maxiter))
    apply_m = (lambda v: v) if precond is None else precond

    nmv = 0
    converged = halted = False
    while nmv < maxiter and not converged and not halted:
        if np.any(x):
            r = apply_m(b - matvec(x))
            nmv += 1
            if callback is not None and callback(float(np.sqrt(r @ r))):
                halted = True
                break
        else:
            r = apply_m(b.copy())
        beta = float(np.sqrt(r @ r))
        if beta <= tol_abs:
            converged = True
            break
        m = min(m_cap, maxiter - nmv)
        if m <= 0:
            break
        q = np.empty((min(m, 64), n))
        q[0] = r / beta
        cols: list[np.ndarray] = []
        cs: list[float] = []
        sn: list[float] = []
        gam = [beta]
        for j in range(m):
            w = apply_m(matvec(q[j]))
            nmv += 1
            qj = q[:j + 1]
            c1 = qj @ w
            w -= qj.T @ c1
            c2 = qj @ w
            w -= qj.T @ c2
            col = np.empty(j + 2)
            col[:j + 1] = c1 + c2
            hn = float(np.sqrt(w @ w))
            col[j + 1] = hn
            for i in range(j):
                t = cs[i] * col[i] + sn[i] * col[i + 1]
                col[i + 1] = -sn[i] * col[i] + cs[i] * col[i + 1]
                col[i] = t
            d = float(np.hypot(col[j], col[j + 1]))
            c_, s_ = (1.0, 0.0) if d == 0.0 else (col[j] / d, col[j + 1] / d)
            cs.append(c_)
            sn.append(s_)
            col[j] = d
            cols.append(col[:j + 1].copy())
            gam.append(-s_ * gam[j])
            gam[j] = c_ * gam[j]
            res = abs(gam[j + 1])
            stop = bool(callback(res)) if callback is not None else False
            if res <= tol_abs:
                converged = True
                break
            if stop:
                halted = True
                break
            if hn == 0.0:
                break
            if j + 1 < m:
                if j + 1 >= q.shape[0]:
                    grown = np.empty((min(m, max(2 * q.shape[0], j + 2)), n))
                    grown[:j + 1] = q[:j + 1]
                    q = grown
                q[j + 1] = w / hn
        k = len(cols)
        if k:
            rmat = np.zeros((k, k))
            for j, col in enumerate(cols):
                rmat[:j + 1, j] = col
            x += q[:k].T @ _back_substitute(rmat, np.asarray(gam[:k]))
    return x, converged, nmv


def krylov_basis_cap(size: int, tmax: int, max_bytes: float) -> int:
    per = max(1.0, 8.0 * float(size))
    return int(max(32, min(int(tmax), int(max_bytes / per))))


# ---------------------------------------------------------------------------
# the driver
# ---------------------------------------------------------------------------
def run_krylov(solver, callback=None, progress_every: int = 0):
    """GMRES on ``(I - T) u = g``; one sweep per product.

    Sweep accounting: 1 for the affine constant ``g``, ``n`` for the Krylov
    products, 1 closing sweep that turns the converged state back into the
    distribution the run reports.  ``iterations`` counts sweeps, so it is
    directly comparable with the source-iteration driver.

    Stopping: GMRES is pushed to ``target`` (``tol * tol_factor``, floored at
    round-off), but once its residual is under ``enough`` (``tol * 1e-2``)
    it gets a bounded extension only -- digits beyond that are cheap but must
    not eat the budget.  **Convergence is decided on the true transport
    residual after the closing sweep**, never on the iterate step.
    """
    c = solver.case
    tol, tmax = c.iteration.tol, c.iteration.tmax
    want_true = c.iteration.true_residual
    ok, why = applicable(solver)
    if not ok:
        raise ValueError(f"Krylov acceleration does not apply: {why}")

    op = FixedPointOperator(solver)
    residuals: list[float] = []
    masses: list[float] = []
    true_res: list[float] = []
    t0 = time.perf_counter()

    def sweep(u):
        out = op.apply(u)
        masses.append(solver.mom.mass)
        residuals.append(np.nan)
        if want_true:
            true_res.append(np.nan)
        return out

    def close(res):
        residuals[-1] = res
        step = len(residuals)
        if progress_every and step % progress_every == 0:
            print(f"  step {step:8d}  residual {res:.6e}  mass {masses[-1]:.10f}",
                  flush=True)
        if callback is not None:
            callback(step, res, solver)

    x = op.state() if solver.restart_info.get("enabled") else np.zeros(op.n)
    g = sweep(np.zeros(op.n))                    # sweep 1: F(0)
    gnorm = float(np.sqrt(g @ g))
    scale = gnorm or 1.0
    close(1.0)

    def matvec(v):
        return v - (sweep(v) - g)

    precond = None
    if c.scheme.krylov_precond and solver.acc is not None:
        precond = MomentPreconditioner(solver, solver.acc)._apply

    enough = max(tol * _ENOUGH_FACTOR, _FLOOR) * scale
    target = max(tol * c.scheme.krylov_tol_factor, _FLOOR) * scale
    cap = krylov_basis_cap(op.n, tmax, c.scheme.krylov_max_bytes)
    reached = -1

    def on_residual(res):
        nonlocal reached
        close(res / scale)
        n = len(residuals)
        if reached < 0 and res <= enough:
            reached = n
        return reached >= 0 and n - reached > max(64, reached // 4)

    converged = False
    final = float("nan")
    for _ in range(c.scheme.krylov_max_rounds):
        budget = tmax - len(residuals) - 1              # keep the closing sweep
        if budget <= 0:
            break
        reached = -1
        x, _, _ = gmres_full(matvec, g, x0=x, tol_abs=target, maxiter=budget,
                             max_basis=cap, precond=precond, callback=on_residual)
        sweep(x)                                        # closing sweep
        final = float(solver.ctx.true_residual(solver.mom, solver.vdf))
        if want_true:
            true_res[-1] = final
        close(final)
        converged = bool(final < tol) if np.isfinite(final) else gnorm == 0.0
        if converged or len(residuals) >= tmax:
            break
        target = max(target * 1.0e-3, _FLOOR * scale)

    wall = time.perf_counter() - t0
    return solver._record(len(residuals), converged, residuals, masses,
                          true_res if want_true else None, wall,
                          scheme="KRYLOV",
                          extra={"krylov_size": op.n, "krylov_basis_cap": cap,
                                 "krylov_precond": precond is not None,
                                 "transport_residual": final})


# ---------------------------------------------------------------------------
# optional preconditioner (experimental; off by default)
# ---------------------------------------------------------------------------
class MomentPreconditioner:
    """Approximate ``(I - T)^-1`` by the macroscopic moment system.

    Helps in the diffusive limit (``tau_R = 1e-4``: 6 343 -> 840 sweeps with
    restarts) and hurts badly in the ballistic limit (``tau_R = 1``: 13 ->
    24 203), because the moment system degenerates there.  No formulation
    valid across the whole Knudsen range has been found; see
    docs/DEFECT_CORRECTION.md.  Kept for experiments, not used by default.
    """

    def __init__(self, solver, acc, sign=-1.0, blocks=(0, 1, 2)):
        self.s = solver
        self.acc = acc
        self.nd = solver.ndof
        self.nt = solver.n_tris
        self.nmom = 3 * self.nd * self.nt
        self.sign = sign
        self.blocks = blocks
        self.calls = 0
        self._zero_vdf = np.zeros_like(solver.vdf)
        self._work = np.zeros((self.nt, 7 * self.nd))

    def _apply(self, r):
        from .acceleration.local_solve import local_solve

        a = self.acc
        nd, nt = self.nd, self.nt
        k = nd * nt
        self._work[:] = 0.0
        parts = (r[:k], r[k:2 * k], r[2 * k:3 * k])
        for slot, blk in enumerate(self.blocks):
            if blk is None:
                continue
            self._work[:, blk * nd:(blk + 1) * nd] = \
                self.sign * parts[slot].reshape(nd, nt).T
        aa_src = np.ascontiguousarray(
            np.einsum("ipq,iq->pi", a.inv_aa_sol, self._work, optimize=True))
        ffa = a.gsolver.rhs(aa_src, self._zero_vdf, self.s.ctx.cxv,
                            self.s.ctx.cyv, self.s.ctx.domega)
        u_trace = a.gsolver.solve(ffa)
        uq = local_solve(aa_src, a.aa_trace, u_trace, self.s.mesh.tri_faces)
        self.calls += 1
        out = np.array(r, dtype=float)
        out[:self.nmom] = np.concatenate([uq[0:nd].ravel(), uq[nd:2 * nd].ravel(),
                                          uq[2 * nd:3 * nd].ravel()])
        return out
