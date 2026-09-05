"""§6.2 stage level: every array against the Fortran dump, per stage.

This is the primary debugging instrument of the port.  A failure here names
the subroutine that diverged; without it, a discrepancy at iteration 400 is
undebuggable.

The dumps are produced by ``tools/build_fortran.sh`` +
``tools/make_golden.py``; when they are absent the whole module skips.
"""
from __future__ import annotations

import pytest

from ..conftest import GOLDEN, has_dump

pytestmark = pytest.mark.fortran

CIS_DUMPS = ["dump_cis_shipped", "dump_cis_deg1", "dump_cis_deg2",
             "dump_cis_np10na20", "dump_cis_tauN1e-2"]
GSIS_DUMPS = ["dump_gsis_shipped", "dump_gsis_deg1"]


def _compare(name, variant=None, rtol=1e-12):
    from compare_stages import compare

    if not has_dump(name):
        pytest.skip(f"no Fortran dump {name}")
    rep = compare(GOLDEN / name / "dump", iters=(1, 2, 3, 10), rtol=rtol,
                  quiet=True, variant=variant)
    bad = [r for r in rep.rows if r[1] != "OK"]
    assert not bad, "stages differing from the Fortran:\n" + "\n".join(
        f"  {n}: {s} {d}" for n, s, d, _ in bad)
    return rep


@pytest.mark.parametrize("name", CIS_DUMPS)
def test_cis_stages_are_bit_exact(name):
    """Gate 3 is stronger than the proposal asks: every CIS array matches the
    Fortran *bit for bit*, not merely to rtol=1e-12."""
    rep = _compare(name, rtol=1e-12)
    not_exact = [r[0] for r in rep.rows if r[3] != "bit-exact" and r[1] == "OK"]
    assert not not_exact, f"no longer bit-exact: {not_exact}"


@pytest.mark.parametrize("name", GSIS_DUMPS)
def test_gsis_stages_variant_a(name):
    """Gate 4 asks rtol=1e-11; the sparse LU reorders the trace solve relative
    to PARDISO, so this one is agreement-to-tolerance rather than bit-exact."""
    _compare(name, variant="A", rtol=1e-11)


def test_gsis_stages_variant_b():
    """§7.2: the TAU_R-rescaled variant, against its own Fortran build."""
    _compare("dump_gsisB_shipped", variant="B", rtol=1e-10)
