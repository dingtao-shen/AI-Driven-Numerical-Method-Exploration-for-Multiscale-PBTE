"""Config parsing, the control.in converter, and the §7.3 restart contract."""
from __future__ import annotations

import numpy as np
import pytest

from pybte import Case
from pybte.config import case_from_control_in, parse_control_in

from ..conftest import CASES, MESH

CONTROL = """&ITERATION
TOL = 1.0d-8,                  !iteration tolerance
TMAX = 8000000,                !Maximum iteration step
/

&GSIS
ACCFLAG = 1,
/

&VELMSH
NPOLE = 20,
NAZIM = 41,
/

&DG
DEG = 3,
/

&FLOW
Cv = 1.d0,
Vg = 1.d0,
TAU_R = 1.d-3,
TAU_N = 1.d5,
TAU_THR = 1.d0,
/

&FILENAME
FNAME_MSH = './A1_Nx11_Ny11.msh',
/

&N_BC
NBC = 4,
/

&BC
BC_NAME = 'SWall', 'NWall', 'EWall', 'WWall',
BC_PHYID = 11, 12, 13, 14,
BC_TYP = 1, 1, 2, 3,
BC_TEMP = 0.d0, 1.d0, 0.d0, 0.d0,
BC_XOFF = 0.d0, 0.d0, 0.d0, 0.d0,
BC_YOFF = 0.d0, 0.d0, 0.d0, 0.d0,
/
"""


def test_parse_control_in(tmp_path):
    p = tmp_path / "control.in"
    p.write_text(CONTROL)
    g = parse_control_in(p)
    assert g["ITERATION"]["TOL"] == 1e-8
    assert g["ITERATION"]["TMAX"] == 8000000
    assert g["GSIS"]["ACCFLAG"] == 1
    assert g["FLOW"]["TAU_N"] == 1e5
    assert g["BC"]["BC_TYP"] == [1, 1, 2, 3]
    assert g["BC"]["BC_NAME"] == ["SWall", "NWall", "EWall", "WWall"]
    assert g["FILENAME"]["FNAME_MSH"] == "./A1_Nx11_Ny11.msh"


def test_case_from_control_in(tmp_path):
    p = tmp_path / "control.in"
    p.write_text(CONTROL)
    c = case_from_control_in(p)
    assert c.iteration.tol == 1e-8
    assert c.scheme.accflag == 1
    assert c.velmesh.nazim == 40, "NAZIM must be forced even, as Fortran does"
    assert c.ndof_tri == 10 and c.ndof_fc == 4 and c.np_tri == 12
    assert c.tau_c == pytest.approx(1.0 / (1.0 / 1e-3 + 1.0 / 1e5))
    assert [b.type for b in c.boundaries] == [
        "thermalising", "thermalising", "nonthermalising", "periodic"]
    assert c.bc_codes() == [1, 1, 2, 3]


def test_yaml_roundtrip(tmp_path):
    c = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    out = tmp_path / "c.yaml"
    c.to_yaml(out)
    c2 = Case.from_yaml(out)
    assert c2.canonical_json() == c.canonical_json()


def test_yaml_unsigned_exponent_is_coerced(tmp_path):
    """YAML 1.1 parses ``1.0e5`` as a *string*; the Fortran habit writes it
    that way, so the loader has to cope rather than blow up later."""
    (tmp_path / "c.yaml").write_text(
        "flow: {cv: 1.0, vg: 1.0, tau_r: 1.0e-3, tau_n: 1.0e5, tau_thr: 1.0}\n"
        "boundaries: [{name: W, phyid: 11, type: thermalising, temp: 0.0}]\n")
    c = Case.from_yaml(tmp_path / "c.yaml")
    assert isinstance(c.flow.tau_n, float) and c.flow.tau_n == 1e5


def test_config_hash_is_sensitive_and_stable(tmp_path):
    c = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    h = c.config_hash
    assert c.config_hash == h
    c.output.dir = "somewhere/else"
    assert c.config_hash == h, "output paths must not change the hash"
    c.flow.tau_r *= 1.0000001
    assert c.config_hash != h, "physics must change the hash"


def test_rejects_bad_values():
    with pytest.raises(ValueError):
        Case.from_dict({"dg": {"deg": 7}})
    with pytest.raises(ValueError):
        Case.from_dict({"scheme": {"acc_variant": "C"}})
    with pytest.raises(ValueError):
        Case.from_dict({"flow": {"tau_r": -1.0}})
    with pytest.raises(ValueError):
        Case.from_dict({"scheme": {"stabilisation": [1.0, 1.0]}})


def test_unknown_boundary_type():
    b = Case.from_dict({"boundaries": [
        {"name": "x", "phyid": 11, "type": "slippery"}]}).boundaries[0]
    with pytest.raises(ValueError, match="unknown boundary type"):
        _ = b.code


# ---------------------------------------------------------------------------
# §7.3 restart contract
# ---------------------------------------------------------------------------
def _tiny_case(tmp_path, **kw):
    c = Case.from_yaml(CASES / "cavity_tauR1e-3_cis.yaml")
    c.velmesh.npole, c.velmesh.nazim = 4, 4
    c.dg.deg = 1
    c.iteration.tmax = 2
    c.basedir = str(tmp_path)
    c.mesh.file = str(MESH)
    c.output.dir = str(tmp_path / "out")
    for k, v in kw.items():
        setattr(c.restart, k, v)
    return c


def test_stale_restart_file_is_an_error(tmp_path):
    """The Fortran silently reads VDF*.out if it happens to exist, changing
    the initial condition and hence the iteration count with no warning."""
    from pybte import Solver

    c = _tiny_case(tmp_path)
    s = Solver(c)
    path = s.save_restart()
    assert path.exists()
    with pytest.raises(FileExistsError, match="restart file"):
        Solver(_tiny_case(tmp_path))


def test_restart_roundtrip(tmp_path):
    from pybte import Solver

    c = _tiny_case(tmp_path)
    s = Solver(c)
    s.run()
    vdf = s.vdf.copy()
    path = s.save_restart()

    c2 = _tiny_case(tmp_path, enabled=True, path=str(path))
    s2 = Solver(c2)
    assert np.array_equal(s2.vdf, vdf)
    assert s2.restart_info["enabled"] is True
    assert len(s2.restart_info["sha256"]) == 64


def test_restart_size_mismatch_is_caught(tmp_path):
    from pybte import Solver

    bad = tmp_path / "junk.out"
    bad.write_bytes(b"\x00" * 64)
    c = _tiny_case(tmp_path, enabled=True, path=str(bad))
    with pytest.raises(ValueError, match="expected"):
        Solver(c)
