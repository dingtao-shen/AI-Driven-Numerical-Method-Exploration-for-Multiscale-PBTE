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
def test_nonthermalising_wall_conserves_energy():
    """A diffusely reflecting wall emits exactly what makes its net normal
    heat flux vanish, so the flux through it must be zero to round-off --
    not merely to discretisation error."""
    s = Solver(_case("cavity_adiabatic_gsis.yaml"))
    s.run()
    q = boundary_heat_flux(s)
    codes = s.bcdata.bc_type[np.clip(s.mesh.face_bc, 0, None)]
    adiabatic = (s.mesh.face_bc >= 0) & (codes == 2)
    thermalising = (s.mesh.face_bc >= 0) & (codes == 1)
    scale = np.abs(q[thermalising]).max()

    assert np.abs(q[adiabatic]).max() < 1e-10 * scale, (
        f"adiabatic faces leak {np.abs(q[adiabatic]).max():.3e} "
        f"against a wall-flux scale of {scale:.3e}")
    # and the two active walls balance: what enters the hot wall leaves the cold
    assert abs(q.sum()) < 1e-8 * scale


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
def test_iterate_residual_understates_the_true_error_for_slow_cis():
    """The pseudo-convergence trap, demonstrated: stop CIS at the default
    tolerance and the true residual is still orders of magnitude larger."""
    c = _case("cavity_tauR1e-3_cis.yaml", **{"flow.tau_r": 1e-2,
                                             "iteration.tol": 1e-6,
                                             "iteration.tmax": 200000})
    s = Solver(c)
    rec = s.run()
    true_res = s.ctx.true_residual(s.mom, s.vdf)
    assert rec.final_residual < 1e-6
    assert true_res > 100 * rec.final_residual, (
        f"iterate residual {rec.final_residual:.2e} vs true {true_res:.2e}")


# ---------------------------------------------------------------------------
# §7.7  sweep cycles
# ---------------------------------------------------------------------------
def test_on_cycle_is_configurable():
    c = _case("cavity_tauR1e-3_cis.yaml", **{"scheme.on_cycle": "break"})
    s = Solver(c)
    assert int(s.order.cycles_broken.sum()) == 0, "the shipped mesh is acyclic"
    with pytest.raises(ValueError, match="on_cycle"):
        Case.from_dict({"scheme": {"on_cycle": "sometimes"}})
