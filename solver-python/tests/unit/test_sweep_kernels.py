"""§Phase-5 gate: the jitted kernels must match the plain reference exactly.

An optimisation that quietly changes the arithmetic -- a different LU, a
reassociated reduction, an inverse instead of a solve -- would otherwise slip
through as "close enough" and only show up as a drifting iteration count
thousands of steps later.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import lu_factor, lu_solve

from pybte import Case, Solver
from pybte._kernels import lu_factor_inplace, lu_solve_inplace
from pybte.sweep_reference import sweep_reference

from ..conftest import CASES


@pytest.mark.parametrize("n", [2, 5, 10, 15])
def test_lu_matches_lapack(n):
    """Our DGETF2/DGETRS look-alikes agree with LAPACK bit for bit."""
    rng = np.random.default_rng(42)
    for _ in range(20):
        a = rng.normal(size=(n, n))
        a += n * np.eye(n)                     # keep it comfortably nonsingular
        b = rng.normal(size=n)

        lu, piv = lu_factor(a.copy())
        want = lu_solve((lu, piv), b.copy())

        aa = a.copy()
        p = np.zeros(n, dtype=np.int32)
        info = lu_factor_inplace(aa, p)
        assert info == 0
        bb = b.copy()
        lu_solve_inplace(aa, p, bb)

        assert np.array_equal(aa, lu), "LU factors differ from LAPACK"
        # scipy already converts LAPACK's 1-based IPIV to 0-based row indices
        assert np.array_equal(p, piv), "pivots differ from LAPACK"
        assert np.array_equal(bb, want), "solution differs from LAPACK"


def test_lu_flags_singular():
    a = np.zeros((4, 4))
    a[0, 0] = 1.0
    p = np.zeros(4, dtype=np.int32)
    assert lu_factor_inplace(a, p) != 0


@pytest.mark.slow
def test_jitted_sweep_matches_reference():
    """One sweep from an identical starting state, both kernels."""
    case = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    case.iteration.tmax = 3
    s = Solver(case)

    # give the state something non-trivial to chew on
    s.step()
    s.step()

    vdf_fast = s.vdf.copy()
    vdf_ref = s.vdf.copy()

    s.ctx.sweep(s.mom, vdf_fast)
    sweep_reference(s.ctx, s.mom, vdf_ref)

    assert np.array_equal(vdf_fast, vdf_ref), (
        "jitted sweep differs from the numpy/LAPACK reference; "
        f"max |diff| = {np.abs(vdf_fast - vdf_ref).max():.3e}")


@pytest.mark.slow
def test_precomputed_matches_onthefly():
    """Hoisting the factorisation out of the iteration must not change a bit."""
    case = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    case.iteration.tmax = 5
    a = Solver(case)
    rec_a = a.run()

    case_b = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    case_b.iteration.tmax = 5
    case_b.performance.precompute_inverse = False
    b = Solver(case_b)
    assert b.ctx.mode == "onthefly"
    rec_b = b.run()

    assert np.array_equal(rec_a.residual_history, rec_b.residual_history)
    assert np.array_equal(rec_a.temp, rec_b.temp)
    assert rec_a.factorisation_count == a.ndir * a.n_tris
    assert rec_b.factorisation_count == 5 * b.ndir * b.n_tris


def test_package_works_without_numba():
    """G2: `pip install pybte` with no optional extras must still solve.

    The kernels fall back to their pure-Python twins, which is far slower but
    must give the *same* answer -- the fallback is a real code path, not a
    stub.  Run in a subprocess because the numba decision is made at import.
    """
    import json
    import os
    import subprocess
    import sys
    from pathlib import Path

    script = """
import json, sys
sys.path.insert(0, %r)
from pybte._numba import HAVE_NUMBA
from pybte import Case, Solver
c = Case.from_yaml(%r)
c.flow.tau_r = 1.0; c.velmesh.npole = 6; c.velmesh.nazim = 8; c.dg.deg = 1
c.output.field = c.output.run_record = c.output.runtime_log = False
r = Solver(c).run()
print(json.dumps({"numba": HAVE_NUMBA, "iters": r.iterations, "mass": r.mass}))
""" % (str(Path(__file__).resolve().parents[2]),
       str(CASES / "cavity_tauR1e-3_cis.yaml"))

    def go(**env_extra):
        env = dict(os.environ, **env_extra)
        out = subprocess.run([sys.executable, "-c", script], check=True,
                             capture_output=True, text=True, env=env)
        return json.loads(out.stdout.strip().splitlines()[-1])

    slow = go(PYBTE_NO_NUMBA="1")
    fast = go()
    assert slow["numba"] is False
    assert slow["iters"] == fast["iters"]
    assert slow["mass"] == fast["mass"], "pure-python path disagrees with the jitted one"


def test_performance_settings_change_nothing_numerically():
    """`performance` selects storage and kernels, never physics.  All three
    settings must give identical answers, or the field is a trap."""
    from pybte._numba import HAVE_NUMBA

    def run(**perf):
        c = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
        c.flow.tau_r = 1.0
        c.velmesh.npole, c.velmesh.nazim, c.dg.deg = 6, 8, 1
        c.iteration.tmax = 12
        c.output.field = c.output.run_record = c.output.runtime_log = False
        for k, v in perf.items():
            setattr(c.performance, k, v)
        s = Solver(c)
        return s, s.run()

    base_s, base = run()
    for perf in ({"precompute_inverse": False},
                 {"threads": 1},
                 {"kernel": "numpy"}):
        s, rec = run(**perf)
        assert np.array_equal(rec.residual_history, base.residual_history), perf
        assert np.array_equal(rec.temp, base.temp), perf

    # ... and the settings actually took effect
    assert run(precompute_inverse=False)[0].ctx.mode == "onthefly"
    assert run(kernel="numpy")[0].ctx.kernel == "numpy"
    if HAVE_NUMBA:
        assert base_s.ctx.kernel == "numba"


def test_performance_rejects_nonsense():
    with pytest.raises(ValueError, match="kernel"):
        Case.from_dict({"performance": {"kernel": "cuda"}})
    with pytest.raises(ValueError, match="threads"):
        Case.from_dict({"performance": {"threads": -4}})
