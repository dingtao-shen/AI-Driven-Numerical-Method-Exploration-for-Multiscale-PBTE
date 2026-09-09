"""Run record and field output -- ``Out_Put_Result.f90``.

Three artefacts come out of a run:

``RunRecord``      the structured result the benchmark verifier reads.  It is
                   complete enough that no downstream tool ever has to parse
                   stdout, and it carries ``sweep_count`` so that grading can
                   be hardware-independent.
``RunTime.txt``    the Fortran's append-only one-line-per-run log.  Format is
                   preserved byte-compatibly -- it is already exactly the
                   benchmark's primary metric.
field dump         the 109x109 sampling grid of ``Out_Put_Conduction``, as
                   ``.npz`` (default) and optionally Tecplot ``.dat``.
"""
from __future__ import annotations

import json
import platform
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["RunRecord", "sampling_grid", "locate_points", "sample_fields",
           "write_tecplot", "append_runtime_log"]


# ---------------------------------------------------------------------------
@dataclass
class RunRecord:
    """The machine-readable run interface.  Do not change field names
    casually -- downstream verifiers read them."""
    iterations: int = 0
    converged: bool = False
    residual_history: np.ndarray = field(default_factory=lambda: np.zeros(0))
    residual_true: np.ndarray | None = None
    mass_history: np.ndarray = field(default_factory=lambda: np.zeros(0))
    temp: np.ndarray = field(default_factory=lambda: np.zeros(0))
    qx: np.ndarray = field(default_factory=lambda: np.zeros(0))
    qy: np.ndarray = field(default_factory=lambda: np.zeros(0))
    temp_dofs: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    qx_dofs: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    qy_dofs: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    sweep_count: int = 0
    factorisation_count: int = 0
    wall_clock: float = 0.0
    setup_wall_clock: float = 0.0
    config_hash: str = ""
    scheme: str = "CIS"
    tol: float = 0.0
    tmax: int = 0
    tau_r: float = 0.0
    tau_n: float = 0.0
    tau_c: float = 0.0
    deg: int = 0
    npole: int = 0
    nazim: int = 0
    n_tris: int = 0
    n_faces: int = 0
    mesh_file: str = ""
    restart: dict = field(default_factory=dict)
    environment: dict = field(default_factory=dict)
    diagnostics: dict = field(default_factory=dict)

    # -- helpers -----------------------------------------------------------
    @property
    def final_residual(self) -> float:
        return float(self.residual_history[-1]) if self.residual_history.size else float("nan")

    @property
    def mass(self) -> float:
        return float(self.temp.sum()) if self.temp.size else float("nan")

    @property
    def contraction(self) -> float:
        """Geometric contraction factor fitted to the tail of the residual
        history, or ``nan`` if the history is too short or not decaying.

        Source iteration converges linearly, so the residual settles into
        ``r_{n+1} = rho r_n``.  ``rho`` approaches 1 as the medium becomes
        optically thick, which is both why CIS crawls and why its residual
        stops being informative about its error.
        """
        h = np.asarray(self.residual_history, dtype=float)
        if h.size < 8:
            return float("nan")
        tail = h[max(1, h.size // 2):]
        tail = tail[np.isfinite(tail) & (tail > 0)]
        if tail.size < 4 or tail[-1] >= tail[0]:
            return float("nan")
        rho = (tail[-1] / tail[0]) ** (1.0 / (tail.size - 1))
        return float(rho) if 0.0 < rho < 1.0 else float("nan")

    @property
    def error_estimate(self) -> float:
        """Estimated distance from the fixed point, ``r * rho/(1 - rho)``.

        The stopping criterion measures the *step* between iterates, not the
        *error*.  For a linearly converging sequence the remaining error is
        the sum of all future steps, ``r rho/(1-rho)``, which for ``rho`` near
        1 is orders of magnitude larger than ``r`` itself.

        This is the number to look at before believing a CIS result: at
        ``tau_R = 1e-4`` on the shipped mesh CIS reports a residual of
        ``2.5e-6`` after 200 000 iterations while its answer is off by 89%.
        """
        rho = self.contraction
        if not np.isfinite(rho):
            return float("nan")
        return float(self.final_residual * rho / (1.0 - rho))

    def to_dict(self, include_fields: bool = True) -> dict:
        d = asdict(self)
        for k, v in list(d.items()):
            if isinstance(v, np.ndarray):
                if not include_fields and v.ndim > 1:
                    d[k] = None
                else:
                    d[k] = v.tolist()
        d["final_residual"] = self.final_residual
        d["mass"] = self.mass
        d["contraction"] = self.contraction
        d["error_estimate"] = self.error_estimate
        return d

    def to_json(self, path, include_fields: bool = True) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(include_fields), indent=2))
        return path

    def to_npz(self, path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {k: v for k, v in asdict(self).items()
                   if isinstance(v, np.ndarray)}
        payload["meta"] = np.array(json.dumps(self.to_dict(include_fields=False)))
        np.savez_compressed(path, **payload)
        return path

    @staticmethod
    def environment_fingerprint() -> dict:
        import numpy as _np
        try:
            import scipy
            scipy_v = scipy.__version__
        except ImportError:      # pragma: no cover
            scipy_v = None
        try:
            import numba
            numba_v = numba.__version__
        except ImportError:      # pragma: no cover
            numba_v = None
        return {
            "python": sys.version.split()[0],
            "numpy": _np.__version__,
            "scipy": scipy_v,
            "numba": numba_v,
            "platform": platform.platform(),
            "machine": platform.machine(),
        }


# ---------------------------------------------------------------------------
def sampling_grid(nx: int = 109, ny: int = 109):
    """``Out_Put_Conduction``'s sampling grid.

    A smootherstep-clustered set of 5 points in each of the first and last
    2% of the interval and 99 uniform points in between.  The clustering
    resolves the kinetic boundary layer, which is why it is reproduced rather
    than replaced by a uniform grid.
    """
    if nx == 109 and ny == 109:
        py = np.zeros(ny)
        for j in range(1, 6):
            a = (j - 1.0) / 8.0
            py[j - 1] = a * a * a * (10.0 - 15.0 * a + 6.0 * a * a) * 0.02
            py[ny - j] = 1.0 - py[j - 1]
        for j in range(6, ny - 4):
            py[j - 1] = (j - 5.0) * (1.0 - 0.02) / 100.0 + 0.01
        return py.copy(), py
    return np.linspace(0.0, 1.0, nx), np.linspace(0.0, 1.0, ny)


def locate_points(mesh, px, py):
    """``(nx, ny)`` array of containing element ids, ``-1`` if outside.

    Reproduces the reference's sub-area test with its ``-1e-12`` slack, and
    its first-match-wins scan over elements in index order.
    """
    nxp, nyp = px.size, py.size
    X, Y = np.meshgrid(px, py, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel()], axis=1)
    out = np.full(pts.shape[0], -1, dtype=np.int32)

    xs = mesh.nodes[mesh.tri_nodes, 0]
    ys = mesh.nodes[mesh.tri_nodes, 1]
    xc = np.concatenate([xs, xs[:, :1]], axis=1)
    yc = np.concatenate([ys, ys[:, :1]], axis=1)

    remaining = np.arange(pts.shape[0])
    for k in range(mesh.n_tris):
        if remaining.size == 0:
            break
        p = pts[remaining]
        inside = np.ones(p.shape[0], dtype=bool)
        for il in range(3):
            area = 0.5 * (p[:, 0] * (yc[k, il] - yc[k, il + 1])
                          - p[:, 1] * (xc[k, il] - xc[k, il + 1])
                          + (xc[k, il] * yc[k, il + 1] - xc[k, il + 1] * yc[k, il]))
            inside &= area >= -1e-12
        hit = remaining[inside]
        out[hit] = k
        remaining = remaining[~inside]
    return out.reshape(nxp, nyp)


def sample_fields(mesh, basis, vdf, cxv, cyv, domega, cv, px, py, triid,
                  ndof_tri=None):
    """Evaluate T, qx, qy on the sampling grid.

    Every point is mapped into its own element's reference coordinates and the
    basis is evaluated there, all at once -- the reference does the same thing
    one point at a time, which on the 109x109 grid is 11 881 separate basis
    evaluations.  The arithmetic is identical; only the loop is gone.
    """
    from .integration import _affine_to_reference

    nxp, nyp = px.size, py.size
    nd = basis.ndof_tri if ndof_tri is None else ndof_tri

    # angular moments of the distribution, per element DOF: (n_tris, ndof)
    m0 = np.einsum("d,dil->il", domega, vdf, optimize=True)
    mx = np.einsum("d,dil->il", cxv * domega, vdf, optimize=True)
    my = np.einsum("d,dil->il", cyv * domega, vdf, optimize=True)

    xs = mesh.nodes[mesh.tri_nodes, 0]
    ys = mesh.nodes[mesh.tri_nodes, 1]
    (bA, cA, dA), (fE, gE, hE) = _affine_to_reference(xs, ys)

    X = np.repeat(px, nyp)
    Y = np.tile(py, nxp)
    k = triid.reshape(-1)
    hit = k >= 0
    kh = k[hit]
    xi = bA[kh] * X[hit] + cA[kh] * Y[hit] + dA[kh]
    eta = fE[kh] * X[hit] + gE[kh] * Y[hit] + hE[kh]

    # (n_hit, ndof) -- basis functions of each point's own element
    phi = np.zeros((kh.size, basis.nodfun_tri.shape[0]))
    for c, (a, b) in enumerate(basis.mono):
        phi += basis.nodfun_tri[:, c] * (xi ** float(a) * eta ** float(b))[:, None]

    def gather(coef):
        out = np.zeros(nxp * nyp)
        out[hit] = np.einsum("pl,pl->p", phi, coef[kh], optimize=True)
        return out.reshape(nxp, nyp)

    pT = gather(m0) / cv
    pqx = gather(mx)
    pqy = gather(my)

    z = np.zeros((nxp, nyp))
    return pT, pqx, pqy, z, z.copy(), z.copy()


def write_tecplot(path, px, py, arrays, names):
    """Point-ordered Tecplot ASCII, matching ``Out_Put_Conduction``'s layout."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    nxp, nyp = px.size, py.size
    with path.open("w") as fh:
        fh.write(' VARIABLES=' + ','.join(f'"{n}"' for n in ["x", "y"] + list(names)) + "\n")
        fh.write(f" ZONE I = {nxp} J = {nyp}\n")
        for j in range(nyp):
            for i in range(nxp):
                vals = [px[i], py[j]] + [a[i, j] for a in arrays]
                fh.write(" " + "".join(f"{v:16.6E}" for v in vals) + "\n")
    return path


def append_runtime_log(path, scheme, tau_r, tau_n, n_tris, npole, nazim,
                       step, seconds):
    """``RunTime.txt``, byte-compatible with the Fortran's
    ``FORMAT(1X,A4,2ES10.2,3I6,I8,F20.1)``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Fortran ES10.2 and Python's 10.2E agree character for character here:
    # both give a leading digit, two decimals and a two-digit signed exponent
    # ("  1.00E-03").  Widths differ only if |exponent| >= 100, which cannot
    # happen for a relaxation time.
    line = (f" {scheme:>4s}{tau_r:10.2E}{tau_n:10.2E}"
            f"{n_tris:6d}{npole:6d}{nazim:6d}{step:8d}{seconds:20.1f}\n")
    with path.open("a") as fh:
        fh.write(line)
    return path
