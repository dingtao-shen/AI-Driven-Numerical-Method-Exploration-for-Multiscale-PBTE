"""§7: every deviation from the Fortran, individually tested.

The proposal lists seven allowed deviations. Each one gets a test that would
fail if the deviation were reverted or mis-implemented; §7.8 (indexing, YAML,
output formats) is covered by the unit tests instead.
"""
from __future__ import annotations

import numpy as np
import pytest

from pybte import Case, Solver
from pybte.bc import boundary_heat_flux

from ..conftest import CASES


def _case(name, **kw):
    c = Case.from_yaml(CASES / name)
    c.output.field = c.output.run_record = c.output.runtime_log = False
    for k, v in kw.items():
        obj, _, attr = k.rpartition(".")
        target = c
        for part in obj.split("."):
            target = getattr(target, part)
        setattr(target, attr, v)
    return c


# ---------------------------------------------------------------------------
# §7.1  boundary-condition dispatch
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_thermalising_is_unchanged_by_the_dispatch():
    """All shipped boundaries are type 1, so restoring the dispatch must be a
    strict superset: the shipped case is unaffected."""
    from ..conftest import GOLDEN, has_golden

    if not has_golden("cis_tauR1e-1_deg3"):
        pytest.skip("no Fortran golden run")
    ref = np.atleast_2d(np.loadtxt(GOLDEN / "cis_tauR1e-1_deg3" / "residual_history.txt"))
    rec = Solver(_case("cavity_tauR1e-3_cis.yaml", **{"flow.tau_r": 1e-1})).run()
    assert np.array_equal(rec.residual_history, ref[:, 1])


@pytest.mark.slow
@pytest.mark.parametrize("accflag", [0, 1])
def test_nonthermalising_wall_conserves_energy(accflag):
    """§7.1: a diffusely reflecting wall emits exactly what makes its own net
    normal heat flux vanish, so that flux must be zero to *round-off* -- the
    emission is constructed from the balance, not approximated.

    Both schemes, because the two failure modes are different: the kinetic
    boundary condition lives in the sweep, the macroscopic trace condition in
    the HDG right-hand side.  GSIS needs the tangential projection of the
    trace flux to be stable at all (docs/FORTRAN_ISSUES.md #9).
    """
    s = Solver(_case("cavity_adiabatic_gsis.yaml", **{"scheme.accflag": accflag}))
    rec = s.run()
    assert rec.converged

    q = boundary_heat_flux(s)
    codes = s.bcdata.bc_type[np.clip(s.mesh.face_bc, 0, None)]
    adiabatic = (s.mesh.face_bc >= 0) & (codes == 2)
    thermalising = (s.mesh.face_bc >= 0) & (codes == 1)
    scale = np.abs(q[thermalising]).max()

    leak = np.abs(q[adiabatic]).max() / scale
    assert leak < 1e-10, f"adiabatic faces leak {leak:.2e} relative"

    # The hot and cold walls balance only to discretisation and to whatever
    # tolerance the run stopped at -- unlike the adiabatic condition, global
    # conservation is not built into the scheme.
    assert abs(q.sum()) / scale < 1e-3

    # T = 0 on one wall, 1 on the opposite, adiabatic sides: the problem is
    # one-dimensional and antisymmetric, so int T dA is exactly 1/2.
    assert rec.mass == pytest.approx(0.5, abs=1e-4)


@pytest.mark.slow
def test_periodic_reproduces_a_one_dimensional_solution():
    """With walls in y and a periodic pair in x the answer cannot depend on x,
    and the cross-plane profile must be the 1-D one."""
    from pybte.io_output import locate_points, sample_fields

    s = Solver(_case("channel_periodic_gsis.yaml"))
    rec = s.run()
    assert rec.converged
    assert int(s.bcdata.periodic_flip.sum()) > 0, "no mirrored pair detected"

    px = np.linspace(0.02, 0.98, 21)
    py = np.linspace(0.05, 0.95, 19)
    tri = locate_points(s.mesh, px, py)
    T, qx, qy, *_ = sample_fields(s.mesh, s.basis, s.vdf, s.ctx.cxv, s.ctx.cyv,
                                  s.ctx.domega, 1.0, px, py, tri)
    assert np.ptp(T, axis=0).max() < 1e-3 * np.ptp(T)
    assert np.abs(qx).max() < 1e-3 * np.abs(qy).mean()
    # int T dA over the unit square with T(0)=0, T(1)=1 is 1/2 by antisymmetry
    assert rec.mass == pytest.approx(0.5, abs=1e-5)


def test_symmetry_bc_is_rejected():
    """BC_TYP=4 reaches the mesh reader but is implemented nowhere; treating
    it as a hot wall silently would be worse than refusing."""
    c = _case("cavity_tauR1e-3_cis.yaml")
    c.boundaries[2].type = "symmetry"
    with pytest.raises(NotImplementedError, match="symmetry"):
        Solver(c)


# ---------------------------------------------------------------------------
# §7.2  acceleration variants
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_variant_b_is_selectable_and_diverges():
    """Variant B is ported faithfully, which means reproducing the fact that
    it diverges -- see docs/FORTRAN_ISSUES.md #6."""
    c = _case("cavity_tauR1e-3_gsis.yaml", **{"scheme.acc_variant": "B",
                                              "iteration.tmax": 40})
    s = Solver(c)
    assert s.acc.variant == "B"
    with np.errstate(all="ignore"):
        rec = s.run()
    assert not rec.converged
    assert not np.isfinite(rec.residual_history[-1]) or rec.residual_history[-1] > 0.9


@pytest.mark.slow
def test_variant_a_and_b_build_different_matrices():
    a = Solver(_case("cavity_tauR1e-3_gsis.yaml", **{"scheme.acc_variant": "A"}))
    b = Solver(_case("cavity_tauR1e-3_gsis.yaml", **{"scheme.acc_variant": "B"}))
    assert a.acc.K.shape == b.acc.K.shape
    assert np.abs((a.acc.K - b.acc.K)).max() > 1e-6


# ---------------------------------------------------------------------------
# §7.4  linear-time assembly
# ---------------------------------------------------------------------------
@pytest.mark.fortran
def test_assembled_matrix_matches_the_fortran_csr():
    """COO assembly must reproduce the Fortran's dense-workspace result --
    same sparsity, same values."""
    import sys
    from pathlib import Path

    from ..conftest import GOLDEN, has_dump

    if not has_dump("dump_gsis_shipped"):
        pytest.skip("no Fortran dump")
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
    from dump_fortran_stages import FortranDump

    d = FortranDump(GOLDEN / "dump_gsis_shipped" / "dump")
    s = Solver(_case("cavity_tauR1e-3_gsis.yaml"))
    K = s.acc.K
    ref = d.csr()
    assert K.shape == ref.shape
    assert K.nnz == ref.nnz
    assert np.abs((K - ref)).max() < 1e-12 * np.abs(ref.data).max()


def test_assembly_is_subquadratic_in_face_count():
    """The reference is O(N_FCS^2); ours must not be.  Compare the ratio of
    assembly times against the ratio of face counts squared."""
    import time

    times, faces = [], []
    for mesh in ("A1_Nx11_Ny11", "A1_Nx21_Ny21"):
        c = _case("cavity_tauR1e-3_gsis.yaml", **{"mesh.file": f"../meshes/{mesh}.msh"})
        t0 = time.perf_counter()
        s = Solver(c)
        times.append(time.perf_counter() - t0)
        faces.append(s.mesh.n_faces)
    ratio = times[1] / times[0]
    quadratic = (faces[1] / faces[0]) ** 2
    assert ratio < 0.5 * quadratic, f"setup scaled {ratio:.1f}x, quadratic is {quadratic:.1f}x"


# ---------------------------------------------------------------------------
# §7.5  stabilisation exposed
# ---------------------------------------------------------------------------
def test_stabilisation_is_configurable():
    base = Solver(_case("cavity_tauR1e-3_gsis.yaml"))
    assert np.allclose(base.acc.st, [1.0, 1.0, 1.0]), "default must match the Fortran"

    tweaked = Solver(_case("cavity_tauR1e-3_gsis.yaml",
                           **{"scheme.stabilisation": [1.0, 2.0, 2.0]}))
    assert np.allclose(tweaked.acc.st, [1.0, 2.0, 2.0])
    assert np.abs(tweaked.acc.K - base.acc.K).max() > 1e-9

    scaled = Solver(_case("cavity_tauR1e-3_gsis.yaml", **{"scheme.scale_by_hmin": True}))
    assert scaled.acc.st[1] == pytest.approx(1.0 / scaled.mesh.hmin)
    assert scaled.acc.st[0] == 1.0


@pytest.mark.slow
def test_stabilisation_does_not_move_the_gsis_answer_much():
    """ST is a numerical-flux parameter, not physics: doubling it changes the
    conditioning and the iteration count but not the converged field by much."""
    a = Solver(_case("cavity_tauR1e-3_gsis.yaml", **{"flow.tau_r": 1e-2})).run()
    b = Solver(_case("cavity_tauR1e-3_gsis.yaml", **{"flow.tau_r": 1e-2,
                                                     "scheme.stabilisation": [2.0, 2.0, 2.0]})).run()
    assert np.abs(a.temp - b.temp).max() / np.abs(a.temp).max() < 5e-3


# ---------------------------------------------------------------------------
# §7.6  true residual alongside the iterate residual
# ---------------------------------------------------------------------------
@pytest.mark.slow
def test_true_residual_is_reported_and_does_not_change_the_stopping_rule():
    plain = Solver(_case("cavity_tauR1e-3_cis.yaml", **{"flow.tau_r": 1e-1})).run()
    withtr = Solver(_case("cavity_tauR1e-3_cis.yaml", **{"flow.tau_r": 1e-1,
                                                         "iteration.true_residual": True})).run()
    assert withtr.iterations == plain.iterations, "stopping rule must be untouched"
    assert np.array_equal(withtr.residual_history, plain.residual_history)
    assert withtr.residual_true is not None
    assert plain.residual_true is None
    assert withtr.residual_true[-1] < 1e-6


@pytest.mark.slow
def test_true_residual_tracks_the_iterate_residual_for_cis():
    """For CIS the two residuals measure the same thing, and this pins that
    down so the distinction is not quietly overclaimed again.

    The sweep leaves ``A f = b(moments_before)`` and the moments are then
    recomputed from that ``f``, so the transport residual evaluated afterwards
    is the change in the moments -- exactly what the iterate residual is.
    Where they *do* differ is GSIS, covered by the next test.
    """
    c = _case("cavity_tauR1e-3_cis.yaml", **{"flow.tau_r": 1e-1,
                                             "iteration.tol": 1e-7,
                                             "iteration.tmax": 200000})
    s = Solver(c)
    rec = s.run()
    true_res = s.ctx.true_residual(s.mom, s.vdf)
    assert 0.1 < true_res / rec.final_residual < 10.0, (
        f"iterate {rec.final_residual:.2e} vs true {true_res:.2e}")


@pytest.mark.slow
def test_true_residual_separates_gsis_from_cis():
    """The measurement behind docs/FORTRAN_ISSUES.md #5: converged to the same
    iterate tolerance, CIS satisfies the discrete transport system and GSIS
    does not."""
    cis = Solver(_case("cavity_tauR1e-3_cis.yaml",
                       **{"flow.tau_r": 1e-1, "iteration.tol": 1e-12,
                          "iteration.tmax": 200000}))
    cis.run()
    gsis = Solver(_case("cavity_tauR1e-3_gsis.yaml",
                        **{"flow.tau_r": 1e-1, "iteration.tol": 1e-12,
                           "iteration.tmax": 200000}))
    gsis.run()

    r_cis = cis.ctx.true_residual(cis.mom, cis.vdf)
    r_gsis = gsis.ctx.true_residual(gsis.mom, gsis.vdf)
    assert r_cis < 1e-10, f"CIS is not at the discrete fixed point: {r_cis:.2e}"
    assert r_gsis > 1e-4, f"GSIS unexpectedly reached it: {r_gsis:.2e}"


@pytest.mark.slow
def test_error_estimate_exposes_pseudo_convergence():
    """§7.6's real content: the stopping criterion measures a step, so the
    error has to be estimated from the observed contraction factor.

    Truncate CIS well short of convergence and check that the estimate lands
    within an order of magnitude of the actual distance from the converged
    answer -- and that the raw residual does not.
    """
    ref = Solver(_case("cavity_tauR1e-3_cis.yaml",
                       **{"flow.tau_r": 1e-2, "iteration.tol": 1e-9,
                          "iteration.tmax": 200000})).run()
    assert ref.converged

    rec = Solver(_case("cavity_tauR1e-3_cis.yaml",
                       **{"flow.tau_r": 1e-2, "iteration.tmax": 2000})).run()
    assert not rec.converged
    actual = float(np.abs(rec.temp - ref.temp).max() / np.abs(ref.temp).max())

    assert 0.9 < rec.contraction < 1.0, f"rho = {rec.contraction}"
    assert actual > 100 * rec.final_residual, (
        "the residual should badly understate the error here: "
        f"residual {rec.final_residual:.2e}, actual error {actual:.2e}")
    ratio = rec.error_estimate / actual
    assert 0.1 < ratio < 30.0, (
        f"error estimate {rec.error_estimate:.2e} vs actual {actual:.2e}")


# ---------------------------------------------------------------------------
# §7.7  sweep cycles
# ---------------------------------------------------------------------------
def test_on_cycle_is_configurable():
    c = _case("cavity_tauR1e-3_cis.yaml", **{"scheme.on_cycle": "break"})
    s = Solver(c)
    assert int(s.order.cycles_broken.sum()) == 0, "the shipped mesh is acyclic"
    with pytest.raises(ValueError, match="on_cycle"):
        Case.from_dict({"scheme": {"on_cycle": "sometimes"}})
