# Proposal 1 — Porting `ACC_2D2V_LinearCallawayModel` from Fortran to Python

**Status:** plan of record for the port phase
**Audience:** Claude Code, executing against this repository
**Source of truth:** `fortran-reference/` (read-only; never edit)
**Deliverable:** `solver-python/`, a self-contained, pip-installable, independently usable solver

---

## 0. Read this first

The Fortran code in `fortran-reference/` is a working research solver for the steady
2D2V gray linear Callaway phonon BTE, discretised with nodal DG in space and discrete
ordinates in angle, offering two iteration schemes selected by a single flag:

- `ACCFLAG = 0` — **CIS** (conventional/source iteration)
- `ACCFLAG = 1` — **GSIS** (general synthetic iterative scheme, HDG macroscopic acceleration)

**The port is a fidelity exercise, not a redesign.** Every numerical choice in the Fortran
is to be reproduced exactly unless this document explicitly says otherwise in §7. If you
find yourself thinking "this would be cleaner if…", stop and check §7. If it is not listed
there, reproduce the Fortran behaviour and open an issue instead.

The reason for this strictness: the Python solver will later become the substrate for a
benchmark whose oracle is the GSIS path and whose baseline is the CIS path. If the port
silently changes convergence behaviour, every downstream task is invalidated and the error
will be discovered months later.

---

## 1. Goals and non-goals

### Goals

- **G1.** A Python package `pybte` that reproduces the Fortran algorithm exactly, to the
  tolerances in §6, for both CIS and GSIS.
- **G2.** Self-contained: no Intel compiler, no MKL, no PARDISO, no licensed dependency.
  `pip install -e .` then `python -m pybte run case.yaml` must work on a clean Linux
  container with only PyPI wheels.
- **G3.** Fast enough that a full CIS run at `TAU_R = 1e-3` on the shipped mesh finishes
  inside the benchmark budget (see §5 for targets).
- **G4.** Structured so that the CIS path and the GSIS path are cleanly separable at module
  boundaries, because task construction will need to strip the latter.
- **G5.** Test-covered at the unit, stage, and integration level, with cross-validation
  against the Fortran.

### Non-goals

- Not adding non-gray / spectral / phonon dispersion. Gray linear Callaway only.
- Not adding time dependence. Steady state only.
- Not adding new acceleration schemes. Port what exists.
- Not building the benchmark. That is Proposal 2.
- Not GPU. CPU only.

---

## 2. What the Fortran actually does — verified inventory

Read this section before touching any code. Line counts are from the shipped archive.

| File | Lines | Role |
|---|---|---|
| `Callaway_2D_2V_DG.f90` | 101 | Main driver; the iteration loop |
| `Global_Parameter.f90` | 136 | Constants; namelist reader for `control.in` |
| `USD_Math.f90` | 477 | Gauss–Legendre, triangle quadrature, factorials, FFT shift (unused) |
| `Spatial_Mesh.f90` | 335 | Gmsh 2.2 ASCII reader; face/neighbour topology; periodic pairing |
| `Velocity_Mesh.f90` | 62 | Angular grid, `CX`, `CY`, `DOMEGA` |
| `Basis_Function.f90` | 191 | Nodal basis on reference triangle and reference face |
| `Integration.f90` | 981 | Precomputed element/face integral tensors |
| `Matrix.f90` | 158 | `Init_Triangle_Order`: per-direction topological sweep ordering |
| `Solvers.f90` | 110 | `DG_Solver_VDF`: the transport sweep |
| `Velocity_Distribution.f90` | 211 | VDF storage, moments, residuals, wall flux, restart I/O |
| `Synthetic_Acceleration.f90` | 686 | GSIS: HDG assembly, HoT source, global/local solve, correction |
| `Synthetic_Acceleration1.f90` | 661 | **Variant** of the above — see §7.2 |
| `Boundary_Conditions.f90` | 48 | BC namelist |
| `Out_Put_Result.f90` | 227 | Field output on a sampling grid; analytic Fourier reference |

### 2.1 The governing system

Steady 2D2V gray linear Callaway model. The unknown `VDF(L, TRID, J1, J2)` is the
perturbed distribution: DG coefficient `L` of element `TRID` at polar index `J1`,
azimuthal index `J2`.

Directions and quadrature (`Velocity_Mesh.f90`):

```
THE   = Gauss-Legendre nodes on [0, pi],   reversed          (NPOLE nodes)
PHI   = Gauss-Legendre on [0,pi] and [pi,2pi], each NAZIM/2, reversed
DOMEGA(J1,J2) = sin(THE(J1)) * WTHE(J1) * WPHI(J2)
CX(J1,J2) = Vg * cos(THE(J1))
CY(J1,J2) = Vg * sin(THE(J1)) * cos(PHI(J2))
```

Note `NAZIM` is forced even: `NAZIM = (NAZIM/2)*2`. Note `sum(DOMEGA) = 4*pi`.

Relaxation times: `TAU_C = 1/(1/TAU_R + 1/TAU_N)`.

Transport equation solved per direction:

```
(1/TAU_C) f  +  c · grad f  =  S
S = Cv*T_s/(4 pi TAU_R)
  + (1/TAU_N) * [ Cv*T_s/(4 pi) + (3/(4 pi Vg^2)) * (qx_s*cx + qy_s*cy) ]
```

Moments (`Calculate_Macro_Properties`):

```
Cv * T_s(L,I) = sum_{J1,J2} VDF(L,I,J1,J2) * DOMEGA
qx_s(L,I)     = sum_{J1,J2} CX * VDF * DOMEGA
qy_s(L,I)     = sum_{J1,J2} CY * VDF * DOMEGA
```

Cell averages `Temp(I)`, `Qx(I)`, `Qy(I)` weight by `INT_NODFUNC_TRI(L,I)`.

### 2.2 The element solve (`Solvers.f90`)

For each direction `(J1,J2)`, elements are visited in `TRI_ORDER(:,J1,J2)`. For each
element a dense `NDOF_TRI x NDOF_TRI` system is built and solved by `DGETRF`/`DGETRS`:

```
A_SOL(M,L) = (1/TAU_C) * INT_NODFUNC_TRI_TRI(M,L,TRID)
           - CX * INT_NODFUNC_TRI_TRI_X(TRID,M,L)
           - CY * INT_NODFUNC_TRI_TRI_Y(TRID,M,L)
           + sum_{IL=1..3} 0.5*(speed + |speed|) * INT_NODFUNC_TRI_FC(TRID,IL,M,L)

speed(IL) = CX * TRIANGLES_INF(TRID, 2*IL) + CY * TRIANGLES_INF(TRID, 2*IL+1)
```

`A_SRC` is the scattering source plus the upwind inflow contributions:

- interior face: `- 0.5*(speed - |speed|) * INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,L) * VDF(L,TRIDext,J1,J2)`
- boundary face: `- 0.5*(speed - |speed|) * FW`, with
  `FW = Cv/(4 pi) * BC_TEMP(BCID) * INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,1)`

**Critical:** because `VDF(:,TRIDext,...)` is read after upwind neighbours have already been
overwritten in this same pass, this is a genuine transport sweep (Gauss–Seidel in the
sweep ordering), *not* a Jacobi update. Reproducing this in-place semantics exactly is
mandatory — a Jacobi port will converge differently and silently invalidate every
iteration count.

### 2.3 The sweep ordering (`Matrix.f90::Init_Triangle_Order`)

Per direction, a topological sort: repeatedly emit any element with no remaining
`-1` (incoming) neighbour flag, then clear the reciprocal flag on its neighbours.
Precomputed once at init for all `NPOLE*NAZIM` directions.

### 2.4 GSIS (`Synthetic_Acceleration.f90`)

Seven macroscopic fields per element DOF, stored in `UQ(7*NDOF_TRI, N_TRIS)`:
`T, qx, qy, sigma_xx, sigma_xy, sigma_yx, sigma_yy`. Three trace fields per face DOF in
`U_TRACE(3*NDOF_FC, N_FCS)`.

Per outer iteration:

1. `Calculate_SRC_ACC_HoTfromDVM` — assemble the higher-order-term source from the
   kinetic solution. `PIxx, PIxy, PIyy` are built from `VDF` contracted with the
   x- and y-derivative tensors and the angular moments
   `CX*(5*CX^2-3)`, `CY*(5*CX^2-1)`, `CX*(5*CY^2-1)`, etc., then scaled by `TAU_C/5`
   and premultiplied by `inv_AA_SOL`.
2. `Global_Problem_Solver_ACC` — solve the sparse global trace system (`KKAcomp`, CSR)
   with PARDISO. **The matrix is assembled once in `Init_Acceleration_New` and does not
   change between iterations**; only the right-hand side does.
3. `Local_Problem_Solver_ACC` — recover `UQ` per element:
   `UQ(:,I) = AA_SRC(:,I) + sum_{IL} AA_TRACE(:,:,IL,I) @ U_TRACE(:, face(I,IL))`
4. `Correct_VDF_Calculate_Macro_Properties` — damped correction:

```
tau_loc = TAU_R / TRIANGLES_Hmin(I)
beta    = min(tau_loc, TAU_THR) / tau_loc
l_T  = (T_ACC  - T_VDF ) * beta
l_qx = (qx_ACC - qx_VDF) * beta
l_qy = (qy_ACC - qy_VDF) * beta
VDF += l_T*Cv/(4 pi) + (CX*l_qx + CY*l_qy) * (TAU_C/TAU_N) * 3/(4 pi Vg^2)
```

Stabilisation `ST(1:3) = 1.0` (with `/Hmin` variants commented out — see §7.5).

### 2.5 Convergence and logging

```
RESIDUAL = sqrt( sum_I (Temp(I)-T_OLD(I))^2 / sum_I Temp(I)^2 )
```

Stop when `RESIDUAL < TOL` or `STEP >= TMAX`.

`RunTime.txt` appends: scheme flag, `TAU_R`, `TAU_N`, `N_TRIS`, `NPOLE`, `NAZIM`, `Step`,
wall-clock. **This is already exactly the benchmark's primary metric.** Preserve it.

### 2.6 The shipped case

`control.in`: `TOL=1e-8`, `TMAX=8e6`, `ACCFLAG=0`, `NPOLE=20`, `NAZIM=40`, `DEG=3`,
`Cv=1`, `Vg=1`, `TAU_R=1e-3`, `TAU_N=1e5`, `TAU_THR=1`.
Mesh `A1_Nx11_Ny11.msh`: 121 nodes, unit square, four physical walls (11=S, 12=N, 13=E,
14=W), all thermalising, `BC_TEMP = 0, 1, 0, 0` — i.e. hot north wall, cold elsewhere.

With `Vg=1`, `L=1`, the Knudsen number is `Kn = TAU_R`, so the shipped case at `1e-3` is
already well into the near-diffusive regime where CIS is expected to crawl.

`Out_Put_Result.f90` contains the analytic Fourier (Laplace) series solution for exactly
this configuration — 200-term sine/sinh series. **This is a genuine diffusion-limit
reference and must be ported**; it is the verification anchor for the asymptotic tasks in
Proposal 2.

---

## 3. Target architecture

```
solver-python/
├── pyproject.toml
├── README.md
├── pybte/
│   ├── __init__.py
│   ├── constants.py           # PI, DBL equivalents, tolerances
│   ├── config.py              # dataclasses + YAML loader; replaces control.in
│   ├── mesh/
│   │   ├── gmsh_reader.py     # Spatial_Mesh.f90 §1  (msh 2.2 ASCII)
│   │   ├── topology.py        # faces, neighbours, BC tags, periodic pairing
│   │   └── geometry.py        # normals, areas, Hmin, reference maps
│   ├── quadrature.py          # USD_Math: gauss_legendre, tri_quadrature
│   ├── basis.py               # Basis_Function.f90
│   ├── integration.py         # Integration.f90 — the precomputed tensors
│   ├── velocity.py            # Velocity_Mesh.f90
│   ├── ordering.py            # Matrix.f90::Init_Triangle_Order
│   ├── sweep.py               # Solvers.f90::DG_Solver_VDF   ← hot path
│   ├── moments.py             # macro properties, residuals
│   ├── bc.py                  # Boundary_Conditions.f90 + wall flux
│   ├── acceleration/
│   │   ├── __init__.py
│   │   ├── assembly.py        # Init_Acceleration_New
│   │   ├── hot_source.py      # Calculate_SRC_ACC_HoTfromDVM
│   │   ├── global_solve.py    # Global_Problem_Solver_ACC (scipy splu)
│   │   ├── local_solve.py     # Local_Problem_Solver_ACC
│   │   └── correction.py      # Correct_VDF_Calculate_Macro_Properties
│   ├── analytic.py            # Fourier series reference
│   ├── io_output.py           # field dump, RunTime.txt, JSON run record
│   └── driver.py              # the outer loop + CLI
├── cases/
│   ├── cavity_tauR1e-3_cis.yaml
│   └── ...
├── meshes/
│   ├── A1_Nx11_Ny11.msh       # copied from Fortran
│   └── *.geo                  # gmsh sources for refinement family
├── tests/
│   ├── unit/
│   ├── stage/                 # cross-validation vs Fortran dumps
│   └── integration/
└── tools/
    ├── build_fortran.sh
    ├── dump_fortran_stages.py
    └── compare_stages.py
```

### 3.1 Config format

Replace the Fortran namelist with YAML. One-to-one field mapping, same names uppercased
to lowercase, plus explicit new fields for things that were implicit:

```yaml
iteration:   {tol: 1.0e-8, tmax: 8000000}
scheme:      {accflag: 0}                 # 0=CIS, 1=GSIS
velmesh:     {npole: 20, nazim: 40}
dg:          {deg: 3}
flow:        {cv: 1.0, vg: 1.0, tau_r: 1.0e-3, tau_n: 1.0e5, tau_thr: 1.0}
mesh:        {file: meshes/A1_Nx11_Ny11.msh}
boundaries:
  - {name: SWall, phyid: 11, type: thermalising, temp: 0.0}
  - {name: NWall, phyid: 12, type: thermalising, temp: 1.0}
  - {name: EWall, phyid: 13, type: thermalising, temp: 0.0}
  - {name: WWall, phyid: 14, type: thermalising, temp: 0.0}
restart:     {enabled: false, path: null}  # see §7.3
output:      {dir: out/, field: true, runtime_log: true, run_record: true}
```

Ship a `control_in_to_yaml.py` converter so existing Fortran cases can be replayed.

---

## 4. Phase plan

Each phase has a hard acceptance gate. Do not start phase N+1 until phase N's gate is
green and its tests are committed.

### Phase 0 — Fortran baseline reproducible (0.5 week)

The port cannot be validated without a runnable reference.

1. Make the code build with **gfortran** instead of ifort. The Makefile already has a
   commented gfortran section. Fix any ifort-isms (`REAL*8`, `INTEGER*8`,
   `MKL_SET_NUM_THREADS`, `OMP_SET_NUM_THREADS`).
2. `MKL_SET_NUM_THREADS`/`MKL_SET_DYNAMIC` must become no-ops or be guarded when not
   linking MKL.
3. **The CIS path needs no PARDISO** — `ACCFLAG=0` never calls `Init_PARDISO`. So CIS can
   be built and run immediately against reference LAPACK/OpenBLAS. Do this first.
4. For GSIS, install MKL from conda-forge or pip (`mkl`, `mkl-devel`) and link PARDISO
   with gfortran. If this proves painful, fall back to §6.4.
5. Pin threading: `OMP_NUM_THREADS=1` for all validation runs. The OpenMP sweep loop is
   over directions and is order-independent, but pin anyway so residual histories are
   bit-reproducible.

**Gate 0:** `./DGACC` runs the shipped `control.in` under gfortran, `OMP_NUM_THREADS=1`,
and produces a residual history. Record it as `fortran-reference/golden/cis_tauR1e-3.log`.
Run twice; the two logs must be byte-identical.

### Phase 1 — Instrument the Fortran for stage dumps (0.5 week)

Add a compile-time `-DDEBUG_DUMP` path (or a `DUMPFLAG` namelist entry) writing raw
binary snapshots after each stage of iterations 1, 2, 3, 10:

- `VDF` after `DG_Solver_VDF`
- `T_s, Qx_s, Qy_s` after `Calculate_Macro_Properties`
- `AA_SRC` after `Calculate_SRC_ACC_HoTfromDVM`
- `U_TRACE` after `Global_Problem_Solver_ACC`
- `UQ` after `Local_Problem_Solver_ACC`
- `VDF, T_s` after `Correct_VDF_Calculate_Macro_Properties`

Also dump the setup-time tensors once: `TRI_ORDER`, `DOMEGA/CX/CY`, all
`INT_NODFUNC_*` arrays, `TRIANGLES_INF`, `TRIANGLES_Hmin`, `FACES_TAG`, `inv_AA_SOL`,
`AA_TRACE`, and the assembled CSR triplet arrays.

**This phase is the single highest-value investment in the whole port.** Without stage
dumps, a discrepancy at iteration 400 is undebuggable; with them, it localises to one
subroutine in minutes.

**Gate 1:** dumps exist for both CIS and GSIS on the shipped case, and
`tools/dump_fortran_stages.py` loads them into numpy arrays with correct shapes and
Fortran ordering.

### Phase 2 — Setup-side port (1.5 weeks)

Port, in this order, validating each against its Phase-1 dump before proceeding:

1. `quadrature.py` — Gauss–Legendre, triangle quadrature tables
2. `velocity.py` — `CX`, `CY`, `DOMEGA`
3. `mesh/` — gmsh reader, topology, geometry
4. `basis.py` — nodal basis
5. `integration.py` — the precomputed tensors
6. `ordering.py` — sweep ordering

Ordering caveats: `TRI_ORDER` is not unique (any valid topological order works), so compare
by *validity* (every element's upwind neighbours precede it) rather than by equality — but
**for the reference implementation, reproduce the Fortran's exact tie-breaking** (ascending
element index scan) so that the Gauss–Seidel sweep produces bitwise-identical results.

**Gate 2:** every setup array matches the Fortran dump to `rtol=1e-13, atol=1e-15`.
`TRI_ORDER` matches exactly.

### Phase 3 — CIS path (1 week)

Port `sweep.py`, `moments.py`, `bc.py`, `driver.py`. Get correctness first with a plain
numpy implementation; optimise in Phase 5.

**Gate 3:**
- single sweep from identical input VDF matches Fortran dump to `rtol=1e-12`
- 50-iteration residual history matches to `rtol=1e-10`
- full run to `TOL=1e-8` produces identical iteration count and a `Temp` field matching to
  `rtol=1e-9`
- reproduced for `TAU_R ∈ {1e-1, 1e-2, 1e-3}` and `DEG ∈ {1,2,3}`

### Phase 4 — GSIS path (2 weeks)

Port `acceleration/`. This is the hardest part. Sub-order:

1. `assembly.py` — the HDG local matrices, `inv_AA_SOL`, `AA_TRACE`, `BA_SOL`, and the
   global sparse matrix. Assemble **directly to COO triplets**, not to the dense row-block
   workspace the Fortran uses (see §7.4).
2. `global_solve.py` — `scipy.sparse.linalg.splu` factorised **once** at init; each
   iteration only calls `.solve()` on a new RHS. This mirrors the PARDISO phase separation
   and is essential to performance.
3. `hot_source.py`, `local_solve.py`, `correction.py`

**Gate 4:**
- assembled global matrix matches the Fortran CSR to `rtol=1e-12` after a common
  permutation
- `AA_SRC`, `U_TRACE`, `UQ`, corrected `VDF` each match the Phase-1 dumps at iterations
  1, 2, 3, 10 to `rtol=1e-11`
- full GSIS run matches Fortran iteration count exactly and converged field to `rtol=1e-9`
- **physics gate:** CIS and GSIS converged `Temp` fields agree to `rtol=1e-8` at every
  `TAU_R` tested. This is the fixed-point-preservation property the benchmark will grade
  on, so it must be verified here, not assumed.

### Phase 5 — Performance (1.5 weeks)

Targets in §5. Work in this order, re-running Gate 3 and Gate 4 after each step:

1. **Precompute and store the factorised element operators.** The Fortran rebuilds and
   LU-factorises `A_SOL` for every element, every direction, *every iteration*, even
   though `A_SOL` depends only on geometry, direction, and `TAU_C` — all iteration-
   invariant. Precompute `inv(A_SOL)` once at setup for all `(element, direction)` pairs
   and reduce the sweep to a matvec. Memory: `N_TRIS * NPOLE * NAZIM * NDOF_TRI^2 * 8` B
   (shipped case, `DEG=3`: ~128 MB — acceptable; add a config flag to fall back to
   on-the-fly factorisation for large meshes).
2. **Numba-jit the sweep kernel** with `prange` over directions, mirroring the Fortran
   OpenMP. Keep the pure-numpy version as `sweep_reference.py` and add a test asserting
   the two agree to `rtol=1e-13`.
3. Vectorise moment computation as a single `einsum`/`tensordot` over `(J1,J2)`.
4. Vectorise the HoT source assembly similarly.

**Gate 5:** performance targets in §5 met; Gates 3 and 4 still green; `sweep.py` and
`sweep_reference.py` agree.

### Phase 6 — Packaging, tests, docs (1 week)

- `pyproject.toml`, `pip install -e .` on a clean container
- CLI: `python -m pybte run case.yaml`, `python -m pybte convert control.in`
- structured run record: JSON with config hash, per-iteration residual, iteration count,
  wall-clock, sweep count, factorisation count, converged moments, environment fingerprint
- `analytic.py` with the Fourier series, plus an L2-error utility against it
- mesh family: `.geo` files and generated `.msh` at 4 refinement levels
- README with equations, discretisation summary, usage, and validation table

**Gate 6:** clean-container install-and-run; `pytest` green; a `VALIDATION.md` table
showing Fortran-vs-Python iteration counts and field errors across the full parameter
sweep in §6.5.

**Total: 8 weeks of focused work.** Phases 3 and 4 are the irreducible core; 0, 1, 2 are
prerequisites that people are tempted to skip and must not.

---

## 5. Performance targets

Measured on the shipped mesh (`N_TRIS` ≈ 200), `DEG=3`, `NPOLE=20`, `NAZIM=40`, single
core unless stated.

| Quantity | Target |
|---|---|
| Setup (mesh, basis, integration, ordering, operator precompute) | < 30 s |
| GSIS setup (HDG assembly + `splu`) | < 30 s |
| Per CIS iteration | < 50 ms |
| Per GSIS iteration | < 150 ms |
| GSIS full run, `TAU_R=1e-3`, `TOL=1e-8` | < 60 s |
| CIS full run, `TAU_R=1e-1`, `TOL=1e-8` | < 5 min |
| Peak RSS | < 2 GB |

Rationale: Proposal 2 requires each benchmark task to complete verification in 10–15
minutes, inside which an agent will run the solver many times. If a single GSIS run costs
more than a minute on the reference case, the task budget cannot be met.

**CIS at `TAU_R=1e-4` is expected to be effectively non-convergent within `TMAX`.** That is
not a bug — it is the phenomenon the whole project is about. Record the iteration count at
which it is truncated.

---

## 6. Validation contract

### 6.1 Unit level

- Gauss–Legendre nodes/weights against `numpy.polynomial.legendre.leggauss` (after
  interval mapping and the Fortran's reversal)
- `sum(DOMEGA) == 4*pi` to `1e-14`
- triangle quadrature integrates degree-`p` polynomials exactly
- nodal basis: partition of unity, Kronecker delta at nodes
- mass matrix row sums equal element area
- mesh: face count identity, every interior face has exactly two neighbours, boundary
  faces carry a valid `BCID`

### 6.2 Stage level

For iterations 1, 2, 3, 10 on the shipped case, compare every array listed in Phase 1
against the Fortran dump. This is the primary debugging instrument. Failures must be
localised to a single stage before any further work.

### 6.3 Integration level

Full runs, comparing iteration count (exact) and converged field (`rtol` per gate).

### 6.4 If PARDISO cannot be linked

Fall back to structural validation of the GSIS path:

- assembled global matrix compared against the Fortran's dumped CSR arrays (this does not
  require *solving* in Fortran, only assembling — so instrument `Init_Acceleration_New`
  and stop there)
- `U_TRACE` validated by residual: `||K @ U_TRACE - F||` computed in Python using the
  Fortran-dumped `K` and `F`
- the fixed-point property (CIS and GSIS agreeing) as the end-to-end check
- published GSIS behaviour as a sanity band: iteration counts should stay roughly flat as
  `TAU_R` decreases, versus CIS growing by orders of magnitude

This fallback is acceptable. Do not let a PARDISO build problem block the port.

### 6.5 Validation sweep

Run the full cross-validation over:

- `TAU_R ∈ {1e0, 1e-1, 1e-2, 1e-3, 1e-4}`
- `TAU_N ∈ {1e5, 1e0, 1e-2}` (RTA-like through to hydrodynamic)
- `DEG ∈ {1, 2, 3}`
- `NPOLE/NAZIM ∈ {(10,20), (20,40)}`
- both schemes

Record every cell in `VALIDATION.md`. Cells where CIS does not converge within `TMAX` are
recorded as such — those are the interesting ones.

---

## 7. Deviations from the Fortran — the complete list

Anything not in this list must be reproduced faithfully.

### 7.1 The boundary-condition dispatch is disabled — restore it

In `Solvers.f90`, the BC branch is entirely commented out:

```fortran
!IF (BC_TYP(BCID).EQ.1) THEN
     FW = Cv/4.d0/PI*BC_TEMP(BCID)*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,1)
     A_SRC = A_SRC - 0.5d0*(speed-ABS(speed))*FW
!END IF
!IF (BC_TYP(BCID).EQ.2) THEN   ... non-thermalising ...
!IF (BC_TYP(BCID).EQ.3) THEN   ... periodic ...
```

The thermalising branch is applied unconditionally to *every* boundary face regardless of
its declared type. The non-thermalising (diffusely reflecting) machinery exists in
`Calculate_FLUX_WALL` but is never wired in; the periodic branch is likewise dead, even
though `Spatial_Mesh.f90` does the periodic face pairing.

**Action:** implement all three branches properly in Python, dispatching on `BC_TYP`, and
verify:
- thermalising reproduces the Fortran exactly (since all shipped BCs are type 1, this is
  a strict superset and Gate 3 is unaffected)
- non-thermalising conserves energy: net normal heat flux through a fully adiabatic
  boundary is zero to `1e-10`
- periodic reproduces a 1D cross-plane analytic result

This is a real gap, not a stylistic one, and Proposal 2 has a task built directly on it.

### 7.2 Two acceleration variants — pick one, keep both

`Synthetic_Acceleration.f90` and `Synthetic_Acceleration1.f90` differ by a **rescaling of
the macroscopic momentum and stress equations by `TAU_R`**:

```
variant A (in Makefile):  AA_SOL(q-block, T-block) = -AA*Cv/3
                          AA_SOL(q-block, q-block) = ST(2)*PP + OO/TAU_R
variant B (not built):    AA_SOL(q-block, T-block) = -AA*Cv/3*TAU_R
                          AA_SOL(q-block, q-block) = OO + ST(2)*PP
```

and correspondingly in `AA_TRACE` and the stress blocks. These are algebraically the same
system; they differ in **conditioning as `TAU_R → 0`**, where variant A's `OO/TAU_R` term
blows up.

**Action:** port variant **A** (the one in the Makefile) as the primary. Port variant B as
an alternative selectable by config (`scheme.acc_variant: A|B`). Then run the §6.5 sweep
on both and record where they diverge in iteration count, accuracy, or global-matrix
condition number. This comparison is scientifically interesting in its own right and
becomes a benchmark task in Proposal 2.

### 7.3 Silent restart file — make it explicit

`Init_Velocity_Distribution_Function` opens `VDF<...>.out` and, if it exists, silently
reads it as the initial condition. For benchmarking this is a determinism hazard: a stale
file left in the working directory changes the initial VDF and therefore the iteration
count, with no warning.

**Action:** restart must be off by default and explicitly configured. When enabled, log the
file path and its hash into the run record. When disabled, error out if a restart file is
present, rather than ignoring it silently.

### 7.4 Global assembly is O(N_FCS²) in time — fix it

`Global_Problem_Solver_ACC` zeroes the dense row-block workspace `KKA(3*NDOF_FC,
N_FCS*3*NDOF_FC)` inside the loop over faces. Memory is fine (O(N_FCS)), but the repeated
zeroing makes assembly quadratic in face count.

**Action:** assemble directly into COO triplet arrays and convert once to CSR. Verify the
resulting matrix equals the Fortran's to `rtol=1e-12`.

### 7.5 Stabilisation parameters — expose them

`ST(1:3) = 1.d0` with `/Hmin` and `/DELTA_REF` scalings commented out. Expose `ST` in the
config (`scheme.stabilisation: [1.0, 1.0, 1.0]`) with an optional `scale_by_hmin` flag, so
the mesh-dependence of the HDG stabilisation can be studied. Default reproduces the
Fortran.

### 7.6 Residual criterion — add a true residual alongside

The Fortran's `RESIDUAL` is the normalised change in cell-average temperature between
successive iterates. For slowly-converging CIS this **systematically underestimates the
true error** (the classic pseudo-convergence trap: a contraction factor near 1 makes
successive iterates close while both are far from the fixed point).

**Action:** keep `Calculate_Residual_T` exactly as-is and report it as `residual_iterate`.
Additionally compute and report `residual_true`: the norm of the actual transport residual
of the discrete system, or in its absence, the error against a converged reference. Both go
in the run record. Do **not** change the stopping criterion — the existing behaviour is the
subject of a benchmark task and must be preserved as the default.

### 7.7 Sweep cycles — handle, do not hang

`Init_Triangle_Order`'s topological loop has no cycle detection. On a non-convex or
sufficiently distorted mesh, some direction can produce a cyclic dependency and the loop
will spin forever emitting nothing.

**Action:** detect a full pass that emits zero elements, and either raise a clear error or
break the cycle by emitting the element with fewest incoming faces (lagging its inflow to
the previous iterate). Log which directions required cycle breaking. Default: raise, so
the shipped meshes are guaranteed acyclic.

### 7.8 Cosmetic, safe to change

- Fortran 1-based indices → Python 0-based, consistently. Document the mapping once in
  `docs/INDEXING.md` and never mix conventions.
- Namelist → YAML (§3.1)
- Tecplot `.dat` output → keep, and additionally emit `.npz` and JSON
- Dead code (`fftshift1d/3d`, `GaussHermit`, `Compress_VDF`) may be dropped

---

## 8. Interfaces the benchmark will depend on

Proposal 2 builds on these; do not change them casually after Phase 6.

```python
from pybte import Case, Solver, RunRecord

case   = Case.from_yaml("cases/cavity_tauR1e-3_cis.yaml")
solver = Solver(case)
record = solver.run()

record.iterations        # int
record.converged         # bool
record.residual_history  # np.ndarray, shape (n_iter,)
record.residual_true     # np.ndarray or None
record.temp              # np.ndarray, shape (N_TRIS,)
record.qx, record.qy
record.temp_dofs         # np.ndarray, shape (NDOF_TRI, N_TRIS)
record.sweep_count       # int  — hardware-independent work metric
record.factorisation_count
record.wall_clock        # float
record.config_hash       # str
record.to_json(path)
```

The run record is what the benchmark verifier reads. It must be complete enough that a
verifier never needs to parse stdout, and it must include `sweep_count` so that grading can
be hardware-independent.

---

## 9. Risk register

| Risk | Mitigation |
|---|---|
| PARDISO unlinkable with gfortran | §6.4 structural fallback; do not block |
| Sweep too slow in Python even with numba | Fall back to a small Cython/C sweep kernel; keep the numpy version as reference. Budget 3 extra days |
| `TRI_ORDER` tie-breaking mismatch causes bitwise divergence | Reproduce the ascending-index scan exactly in Phase 2; assert equality at Gate 2 |
| `inv(A_SOL)` precompute exhausts memory on refined meshes | Config flag to fall back to on-the-fly factorisation; document the crossover |
| Fortran itself has a latent bug we faithfully reproduce | §6.1 unit tests, energy-balance checks, and comparison against the analytic Fourier solution are independent of the Fortran and will surface it. If found, document in `docs/FORTRAN_ISSUES.md` and decide explicitly whether to reproduce or fix |
| Port drifts into redesign | This document, §7 as the closed list of allowed deviations |

---

## 10. Definition of done

- [ ] `pip install -e .` on a clean `python:3.11-slim` container, then a CIS and a GSIS run
      both succeed with no system dependencies beyond PyPI wheels
- [ ] `pytest` green: unit, stage, integration
- [ ] `VALIDATION.md` covers the full §6.5 sweep with Fortran-vs-Python iteration counts
      and field errors
- [ ] CIS and GSIS converge to the same discrete fixed point to `rtol=1e-8` at every
      tested `TAU_R` — **this is the load-bearing property for the whole benchmark**
- [ ] Performance targets in §5 met
- [ ] All §7 deviations implemented and individually tested
- [ ] `RunRecord` interface (§8) stable and documented
- [ ] `docs/` contains the governing equations, the discretisation, the indexing
      convention, and the known-issues list
