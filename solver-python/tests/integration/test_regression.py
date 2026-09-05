"""Frozen results, so the solver cannot drift.

Every number here was verified against the Fortran reference before that
reference was removed from the tree (see VALIDATION.md and the git history at
the `Port ACC_2D2V_LinearCallawayModel` commit).  The iteration counts were
*bit-identical* to the Fortran for CIS and exact for GSIS; the masses are this
solver's own converged values.

These are the guard the Fortran comparison used to be.  A change that alters
any of them is a change in the numerics, and has to be justified rather than
accepted.
"""
from __future__ import annotations

import numpy as np
import pytest

from pybte import Case, Solver

from ..conftest import CASES

pytestmark = pytest.mark.slow

#: (tau_r, tau_n, deg, npole, nazim, scheme) -> (iterations, mass, converged)
#:
#: CIS entries below tau_R = 1e-2 are truncated at 200 000 iterations on
#: purpose: CIS does not converge there, and the mass shows how far off it is
#: (0.0282 against 0.25 at tau_R = 1e-4).
FROZEN = {
    (1.0, 1e5, 3, 20, 40, "CIS"): (24, 0.24968704672388894, True),
    (1e-1, 1e5, 3, 20, 40, "CIS"): (317, 0.24993014541986364, True),
    (1e-2, 1e5, 1, 20, 40, "CIS"): (16598, 0.2499894931642579, True),
    (1e-2, 1e5, 2, 20, 40, "CIS"): (16832, 0.24998995914186684, True),
    (1e-2, 1e5, 3, 10, 20, "CIS"): (16833, 0.2499661216015033, True),
    (1e-2, 1e5, 3, 20, 40, "CIS"): (16836, 0.24998994666784718, True),
    (1e-2, 1e0, 3, 20, 40, "CIS"): (16988, 0.24998991833354237, True),
    (1e-2, 1e-2, 3, 20, 40, "CIS"): (31482, 0.24998679105767513, True),
    (1.0, 1e5, 3, 20, 40, "GSIS"): (22, 0.2496873776739861, True),
    (1e-1, 1e5, 3, 20, 40, "GSIS"): (30, 0.24999647389562388, True),
    (1e-2, 1e5, 1, 20, 40, "GSIS"): (35, 0.24995968584823433, True),
    (1e-2, 1e5, 2, 20, 40, "GSIS"): (29, 0.24995682731416324, True),
    (1e-2, 1e5, 3, 10, 20, "GSIS"): (26, 0.24996286908141993, True),
    (1e-2, 1e5, 3, 20, 40, "GSIS"): (27, 0.24996193435163644, True),
    (1e-2, 1e0, 3, 20, 40, "GSIS"): (27, 0.2499617561539731, True),
    (1e-2, 1e-2, 3, 20, 40, "GSIS"): (31, 0.24995748636168594, True),
    (1e-3, 1e5, 3, 20, 40, "GSIS"): (50, 0.24999752787448726, True),
    (1e-4, 1e5, 3, 20, 40, "GSIS"): (281, 0.24999978596626848, True),
}

#: the cheap subset, for a run that finishes in about a minute
FAST = [k for k in FROZEN if k[0] >= 1e-1 or k[5] == "GSIS"]


def _run(key, tmax=200_000):
    tau_r, tau_n, deg, npole, nazim, scheme = key
    c = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    c.scheme.accflag = 1 if scheme == "GSIS" else 0
    c.flow.tau_r, c.flow.tau_n = tau_r, tau_n
    c.dg.deg = deg
    c.velmesh.npole, c.velmesh.nazim = npole, nazim
    c.iteration.tmax = tmax
    c.output.field = c.output.run_record = c.output.runtime_log = False
    return Solver(c).run()


@pytest.mark.parametrize("key", FAST, ids=lambda k: f"{k[5]}_tau{k[0]:g}_deg{k[2]}_{k[3]}x{k[4]}")
def test_frozen_results(key):
    want_iters, want_mass, want_conv = FROZEN[key]
    rec = _run(key)
    assert rec.converged is want_conv
    assert rec.iterations == want_iters, (
        f"iteration count moved: {rec.iterations} vs {want_iters}")
    assert rec.mass == pytest.approx(want_mass, rel=0, abs=1e-12), (
        f"converged mass moved: {rec.mass!r} vs {want_mass!r}")


@pytest.mark.parametrize("key", [k for k in FROZEN if k not in FAST],
                         ids=lambda k: f"{k[5]}_tau{k[0]:g}_deg{k[2]}_{k[3]}x{k[4]}")
def test_frozen_results_slow_cis(key):
    """The CIS cells that take minutes each -- same assertion, split out so
    the fast set can be run alone."""
    test_frozen_results(key)


def test_cis_iteration_count_grows_as_kn_falls():
    """The phenomenon the project exists to measure, frozen as a shape rather
    than as individual numbers: 24 -> 317 -> 16836 for CIS while GSIS stays
    at 22 -> 30 -> 27."""
    cis = [FROZEN[(t, 1e5, 3, 20, 40, "CIS")][0] for t in (1.0, 1e-1, 1e-2)]
    gsis = [FROZEN[(t, 1e5, 3, 20, 40, "GSIS")][0] for t in (1.0, 1e-1, 1e-2)]
    assert cis[1] > 10 * cis[0]
    assert cis[2] > 50 * cis[1]
    assert max(gsis) < 2 * min(gsis)
    assert cis[2] / gsis[2] > 500


def test_field_matches_the_frozen_sample():
    """A few field values on the 109x109 output grid, frozen from the run that
    matched the Fortran's own Tecplot output to 5e-8 (its print precision)."""
    from pybte.io_output import locate_points, sample_fields, sampling_grid

    c = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    c.flow.tau_r = 1e-1
    c.output.field = c.output.run_record = c.output.runtime_log = False
    s = Solver(c)
    s.run()
    px, py = sampling_grid()
    tri = locate_points(s.mesh, px, py)
    assert not np.any(tri < 0)
    T, qx, qy, *_ = sample_fields(s.mesh, s.basis, s.vdf, s.ctx.cxv, s.ctx.cyv,
                                  s.ctx.domega, c.flow.cv, px, py, tri)
    # The grid is clustered, so px[0] is just inside the cold corner rather
    # than on it -- hence 3.3e-3 and not 0.
    assert T[0, 0] == pytest.approx(0.0032735050690818894, rel=1e-9)
    assert T[54, 54] == pytest.approx(0.2499447278527214, rel=1e-9)
    assert T[54, 104] == pytest.approx(0.8728271107454624, rel=1e-9)
    assert float(T.min()) == pytest.approx(0.0032735050690818894, rel=1e-9)
    assert float(T.max()) == pytest.approx(0.8983932311211115, rel=1e-9)
    assert float(np.abs(qx).max()) == pytest.approx(0.1883719557037961, rel=1e-9)
    assert float(T.mean()) == pytest.approx(0.24989294601937023, rel=1e-9)
