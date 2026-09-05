#!/usr/bin/env bash
# build_fortran.sh -------------------------------------------------------
# Build the Fortran reference solver with gfortran.
#
# The tree under fortran-reference/ is the source of truth and is never
# modified.  We copy it into fortran-build/src/<flavour>/, overlay the
# additive shims from tools/fortran_shims/, and compile there.
#
# Usage:
#   tools/build_fortran.sh [--variant A|B] [--no-mkl] [--plain] [--openmp]
#
#   --variant A   Synthetic_Acceleration.f90   (default, the Makefile's choice)
#   --variant B   Synthetic_Acceleration1.f90  (the TAU_R-rescaled variant)
#   --no-mkl      link reference LAPACK/BLAS and use the dense PARDISO stub
#   --plain       build the *unmodified* reference driver (no dumps, and
#                 OMP_SET_NUM_THREADS(24) as shipped)
#   --openmp      compile with -fopenmp (NOT recommended: the reference
#                 privatises module ALLOCATABLEs, and we want bit-repro)
#
# Output: fortran-build/<flavour>/DGACC
#-------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
REF="$ROOT/fortran-reference/ACC_2D2V_LinearCallawayModel"
SHIM="$HERE/fortran_shims"
BUILD="$ROOT/fortran-build"

VARIANT=A
USE_MKL=1
PLAIN=0
OPENMP=0
NO_FIXES=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --variant)  VARIANT="$2"; shift 2 ;;
    --no-mkl)   USE_MKL=0; shift ;;
    --plain)    PLAIN=1; shift ;;
    --openmp)   OPENMP=1; shift ;;
    --no-fixes) NO_FIXES=1; shift ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done

FLAVOUR="var${VARIANT}"
[[ $USE_MKL  -eq 0 ]] && FLAVOUR="${FLAVOUR}_nomkl"
[[ $PLAIN   -eq 1 ]] && FLAVOUR="${FLAVOUR}_plain"
[[ $OPENMP  -eq 1 ]] && FLAVOUR="${FLAVOUR}_omp"
[[ $NO_FIXES -eq 1 ]] && FLAVOUR="${FLAVOUR}_nofix"

SRC="$BUILD/$FLAVOUR/src"
BIN="$BUILD/$FLAVOUR/bin"
rm -rf "$BUILD/$FLAVOUR"
mkdir -p "$SRC" "$BIN"

# ---- assemble the source set -------------------------------------------
for f in Global_Parameter.f90 USD_Math.f90 Boundary_Conditions.f90 \
         Spatial_Mesh.f90 Velocity_Mesh.f90 Matrix.f90 Basis_Function.f90 \
         Integration.f90 Velocity_Distribution.f90 Solvers.f90 \
         Out_Put_Result.f90 Callaway_2D_2V_DG.f90; do
  cp "$REF/$f" "$SRC/$f"
done

case "$VARIANT" in
  A) cp "$REF/Synthetic_Acceleration.f90"  "$SRC/Synthetic_Acceleration.f90" ;;
  B) cp "$REF/Synthetic_Acceleration1.f90" "$SRC/Synthetic_Acceleration.f90" ;;
  *) echo "variant must be A or B" >&2; exit 1 ;;
esac

# ---- correctness fixes applied to the build copy ------------------------
# FIX-1 (docs/FORTRAN_ISSUES.md #1): Init_Acceleration_New allocates
#   AA_TMP(3*NDOF_FC,3*NDOF_FC) and only ever assigns its three diagonal
#   NDOF_FC-blocks.  The six off-diagonal blocks are never initialised, so
#   they contain whatever the allocator handed back -- in practice the
#   just-freed AA_SOL workspace, entries of order 1e5.  Those leak straight
#   into every row of the HDG global matrix and destroy the GSIS solution.
#   This is undefined behaviour, not a numerical choice, so there is nothing
#   to "port faithfully"; we zero it.
if [[ $NO_FIXES -eq 0 ]]; then
  python3 - "$SRC/Synthetic_Acceleration.f90" <<'PY'
import sys, pathlib
p = pathlib.Path(sys.argv[1]); s = p.read_text()
anchor = "ALLOCATE(C1A_TMP(3*NDOF_FC),C2A_TMP(3*NDOF_FC))"
assert s.count(anchor) == 1, "anchor for FIX-1 not found exactly once"
s = s.replace(anchor, anchor + "\nAA_TMP = 0.d0   !! FIX-1: see tools/build_fortran.sh")
p.write_text(s)
print("  FIX  AA_TMP zero-initialised in Synthetic_Acceleration.f90")
PY
fi

cp "$SHIM/omp_stub.f90" "$SRC/"
OBJ_EXTRA=""

if [[ $PLAIN -eq 0 ]]; then
  cp "$SHIM/Dump_Utils.f90" "$SRC/"
  cp "$SHIM/Callaway_2D_2V_DG_instrumented.f90" "$SRC/Callaway_2D_2V_DG.f90"
  OBJ_EXTRA="Dump_Utils"
fi

if [[ $USE_MKL -eq 0 ]]; then
  cp "$SHIM/mkl_stub.f90" "$SRC/"
  cp "$SHIM/pardiso_stub.f90" "$SRC/"
fi

# ---- compiler flags -----------------------------------------------------
FC=${FC:-gfortran}
# -O2 without -ffast-math: no FP reassociation, so the residual history is
# bit-reproducible run to run.  -fallow-argument-mismatch is needed because
# Spatial_Mesh.f90 reads into differently-shaped slices of the same array.
FFLAGS="-O2 -J$BIN -I$BIN -fallow-argument-mismatch -std=legacy -w"
[[ $OPENMP -eq 1 ]] && FFLAGS="$FFLAGS -fopenmp"

MKLROOT_GUESS="${CONDA_PREFIX:-/home/dtshen/miniconda3/envs/dev}"
if [[ $USE_MKL -eq 1 ]]; then
  MKLLIB="${MKL_LIB_DIR:-$MKLROOT_GUESS/lib}"
  if [[ ! -f "$MKLLIB/libmkl_gf_lp64.so" ]]; then
    echo "MKL not found in $MKLLIB -- re-run with --no-mkl" >&2
    exit 1
  fi
  LDFLAGS="-L$MKLLIB -Wl,-rpath,$MKLLIB -Wl,--no-as-needed \
           -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
else
  LDFLAGS="-llapack -lblas"
fi
[[ $OPENMP -eq 1 ]] && LDFLAGS="$LDFLAGS -fopenmp"

# ---- compile ------------------------------------------------------------
ORDER=(omp_stub Global_Parameter USD_Math Boundary_Conditions Spatial_Mesh \
       Velocity_Mesh Matrix Basis_Function Integration Velocity_Distribution \
       Solvers Synthetic_Acceleration Out_Put_Result)
[[ -n "$OBJ_EXTRA" ]] && ORDER+=("$OBJ_EXTRA")
[[ $USE_MKL -eq 0 ]] && ORDER+=(mkl_stub pardiso_stub)
ORDER+=(Callaway_2D_2V_DG)

OBJS=()
for m in "${ORDER[@]}"; do
  echo "  FC  $m.f90"
  $FC $FFLAGS -c "$SRC/$m.f90" -o "$BIN/$m.o"
  OBJS+=("$BIN/$m.o")
done

echo "  LD  DGACC"
$FC "${OBJS[@]}" $LDFLAGS -o "$BUILD/$FLAVOUR/DGACC"
echo "built: $BUILD/$FLAVOUR/DGACC"
