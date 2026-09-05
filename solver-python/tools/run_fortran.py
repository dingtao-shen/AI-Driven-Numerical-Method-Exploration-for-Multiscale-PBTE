#!/usr/bin/env python3
"""Run the Fortran reference solver on a generated ``control.in``.

The reference reads ``./control.in`` and ``FNAME_MSH`` relative to the
current directory, so every run gets its own scratch directory.  This
script is the single entry point used by

  * the Phase-0 golden runs,
  * the Phase-1 stage dumps,
  * the §6.5 validation sweep.

Usage
-----
    python tools/run_fortran.py --out DIR [--tau-r 1e-3] [--accflag 0] ...

Writes into DIR:  control.in, the mesh, stdout.log, residual_history.txt,
RunTime.txt, the Tecplot field files, and (with --dump) the raw stage
dumps in DIR/dump/.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "fortran-reference" / "ACC_2D2V_LinearCallawayModel"
BUILD = ROOT / "fortran-build"

CONTROL_TEMPLATE = """&ITERATION
TOL = {tol},                  !iteration tolerance
TMAX = {tmax},                !Maximum iteration step
/

&GSIS
ACCFLAG = {accflag},               !0-CIS, 1-GSIS
/

&VELMSH
NPOLE = {npole},                      !# of discrete polar angle
NAZIM = {nazim},                      !# of discrete azimuthal angle
/

&DG
DEG = {deg},                       ! degree of polynomial in DG discretization
/

&FLOW
Cv = {cv},                    !specific heat
Vg = {vg},                    !magnitude of group velocity
TAU_R = {tau_r},                      !mean free time of resistive scattering
TAU_N = {tau_n},                     !mean free time of normal scattering
TAU_THR = {tau_thr},                    !threhold value in GSIS correction
/

&FILENAME
FNAME_MSH = './{mesh}',      !file name for spatial mesh
/

&N_BC
NBC = {nbc},
/

&BC
BC_NAME = {bc_name},
BC_PHYID = {bc_phyid},
BC_TYP = {bc_typ},
BC_TEMP = {bc_temp},
BC_XOFF = {bc_xoff},
BC_YOFF = {bc_yoff},
/
"""


def fortran_double(x: float) -> str:
    """1e-3 -> 1.0d-3 ; the reference namelist reader wants d-exponents."""
    return repr(float(x)).replace("e", "d") if "e" in repr(float(x)) else f"{float(x)!r}d0"


def build_control(args) -> str:
    n = len(args.bc_phyid)
    return CONTROL_TEMPLATE.format(
        tol=fortran_double(args.tol),
        tmax=args.tmax,
        accflag=args.accflag,
        npole=args.npole,
        nazim=args.nazim,
        deg=args.deg,
        cv=fortran_double(args.cv),
        vg=fortran_double(args.vg),
        tau_r=fortran_double(args.tau_r),
        tau_n=fortran_double(args.tau_n),
        tau_thr=fortran_double(args.tau_thr),
        mesh=Path(args.mesh).name,
        nbc=n,
        bc_name=", ".join(f"'{s}'" for s in args.bc_name),
        bc_phyid=", ".join(str(v) for v in args.bc_phyid),
        bc_typ=", ".join(str(v) for v in args.bc_typ),
        bc_temp=", ".join(fortran_double(v) for v in args.bc_temp),
        bc_xoff=", ".join(fortran_double(v) for v in args.bc_xoff),
        bc_yoff=", ".join(fortran_double(v) for v in args.bc_yoff),
    )


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", required=True)
    p.add_argument("--exe", default=None, help="path to DGACC (default fortran-build/varA/DGACC)")
    p.add_argument("--variant", default="A", choices=["A", "B"])
    p.add_argument("--mesh", default=str(REF / "A1_Nx11_Ny11.msh"))
    p.add_argument("--tol", type=float, default=1e-8)
    p.add_argument("--tmax", type=int, default=8_000_000)
    p.add_argument("--accflag", type=int, default=0)
    p.add_argument("--npole", type=int, default=20)
    p.add_argument("--nazim", type=int, default=40)
    p.add_argument("--deg", type=int, default=3)
    p.add_argument("--cv", type=float, default=1.0)
    p.add_argument("--vg", type=float, default=1.0)
    p.add_argument("--tau-r", type=float, default=1e-3)
    p.add_argument("--tau-n", type=float, default=1e5)
    p.add_argument("--tau-thr", type=float, default=1.0)
    p.add_argument("--bc-name", nargs="+", default=["SWall", "NWall", "EWall", "WWall"])
    p.add_argument("--bc-phyid", nargs="+", type=int, default=[11, 12, 13, 14])
    p.add_argument("--bc-typ", nargs="+", type=int, default=[1, 1, 1, 1])
    p.add_argument("--bc-temp", nargs="+", type=float, default=[0.0, 1.0, 0.0, 0.0])
    p.add_argument("--bc-xoff", nargs="+", type=float, default=[0.0, 0.0, 0.0, 0.0])
    p.add_argument("--bc-yoff", nargs="+", type=float, default=[0.0, 0.0, 0.0, 0.0])
    p.add_argument("--dump", action="store_true")
    p.add_argument("--dump-iters", default=None)
    p.add_argument("--dump-no-vdf", action="store_true")
    p.add_argument("--timeout", type=float, default=None)
    args = p.parse_args(argv)

    exe = Path(args.exe) if args.exe else BUILD / f"var{args.variant}" / "DGACC"
    if not exe.exists():
        print(f"missing executable {exe}; run tools/build_fortran.sh first", file=sys.stderr)
        return 2

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    shutil.copy(args.mesh, out / Path(args.mesh).name)
    (out / "control.in").write_text(build_control(args))

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    if args.dump:
        dumpdir = (out / "dump").resolve()
        dumpdir.mkdir()
        env["PYBTE_DUMP_DIR"] = str(dumpdir)
        if args.dump_iters:
            env["PYBTE_DUMP_ITERS"] = args.dump_iters
        if args.dump_no_vdf:
            env["PYBTE_DUMP_NO_VDF"] = "1"
    else:
        env.pop("PYBTE_DUMP_DIR", None)

    t0 = time.time()
    with open(out / "stdout.log", "w") as fh:
        proc = subprocess.run([str(exe.resolve())], cwd=out, env=env,
                              stdout=fh, stderr=subprocess.STDOUT,
                              timeout=args.timeout)
    dt = time.time() - t0
    (out / "wallclock.txt").write_text(f"{dt:.3f}\n")
    print(f"[run_fortran] {out}  rc={proc.returncode}  {dt:.1f}s")
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
