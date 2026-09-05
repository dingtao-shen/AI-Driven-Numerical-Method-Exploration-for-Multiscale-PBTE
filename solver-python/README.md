# pybte

Steady 2D2V gray linear **Callaway phonon Boltzmann** solver: nodal
discontinuous Galerkin in space, discrete ordinates in angle, with two
iteration schemes —

* **CIS** — conventional source iteration (`scheme.accflag: 0`)
* **GSIS** — general synthetic iterative scheme, HDG macroscopic acceleration
  (`scheme.accflag: 1`)

This is a fidelity port of the Fortran research solver in
`../fortran-reference/ACC_2D2V_LinearCallawayModel`. Self-contained: no Intel
compiler, no MKL, no PARDISO — `pip install -e .` and PyPI wheels only.

## Install and run

```bash
pip install -e ".[dev]"          # numba + pytest; plain install works without
python -m pybte info cases/cavity_tauR1e-3_cis.yaml
python -m pybte run  cases/cavity_tauR1e-3_gsis.yaml
python -m pybte convert ../fortran-reference/ACC_2D2V_LinearCallawayModel/control.in
```

`numba` is optional. Without it every kernel falls back to its pure-Python
twin — correct but far slower; the test suite checks the two agree.

## API

```python
from pybte import Case, Solver

case   = Case.from_yaml("cases/cavity_tauR1e-3_gsis.yaml")
solver = Solver(case)
record = solver.run()

record.iterations          # int
record.converged           # bool
record.residual_history    # (n_iter,)  the Fortran's iterate residual
record.residual_true       # (n_iter,) or None -- a real residual, see §7.6
record.temp, record.qx, record.qy        # (n_tris,) element integrals
record.temp_dofs                          # (ndof_tri, n_tris)
record.sweep_count         # hardware-independent work metric
record.factorisation_count
record.wall_clock, record.config_hash
record.to_json("run.json")
```

`Temp`, `Qx` and `Qy` are element **integrals**, not averages: `record.temp.sum()`
is the domain integral of `T` (0.25 for the shipped cavity). Divide by
`solver.mesh.tri_area` for cell averages.

## The physics in one paragraph

Phonons stream and scatter; Callaway splits scattering into a resistive
channel (`tau_R`, destroys momentum) and a normal one (`tau_N`, conserves it).
With `Vg = 1` on a unit domain the Knudsen number is just `tau_R`. As
`tau_R → 0` the model reduces to Fourier conduction with
`kappa = Cv Vg² tau_R / 3`, and the shipped cavity's diffusion limit is the
Laplace solution carried in `pybte.analytic`. See `docs/EQUATIONS.md` for the
discretisation.

## Why GSIS exists

Source iteration regenerates almost all of its own source when the medium is
optically thick, so its spectral radius approaches 1 as `tau_R → 0`. Measured
on the shipped 200-element mesh at `DEG=3`:

| `tau_R` | CIS iterations | GSIS iterations |
|---|---|---|
| 1e+00 | 24 | 22 |
| 1e-01 | 317 | 30 |
| 1e-02 | 16 836 | 27 |
| 1e-03 | > 20 000 (not converged) | 50 |

**CIS at `tau_R = 1e-4` is expected to be effectively non-convergent within any
practical budget.** That is not a bug; it is the phenomenon the project is
about.

## Fidelity

Verified against the Fortran on every stage of every iteration
(`tests/stage/`, `tools/compare_stages.py`):

* **Setup — bit-identical.** Mesh topology, geometry, angular quadrature,
  both nodal bases, all nine integral tensors and `TRI_ORDER` reproduce the
  Fortran exactly, not merely to tolerance. That required matching the
  reference's summation orders and its use of libm `pow` for monomials.
* **CIS — bit-identical over full runs**, the longest verified being 16 836
  iterations at `tau_R = 1e-2`: identical iteration count, and every entry of
  the residual and mass histories equal to the last bit. The element LU is a
  hand-written `DGETF2`/`DGETRS` look-alike for exactly this reason.
* **GSIS — agrees to ~2.5e-13** relative, stage by stage. The residual is the
  sparse LU reordering the trace solve relative to PARDISO. The assembled
  global matrix matches the Fortran CSR to 1.5e-15 with identical sparsity.

`VALIDATION.md` has the full §6.5 sweep.

## Two things you should know before trusting a result

Both are documented in full in **`docs/FORTRAN_ISSUES.md`**.

**1. The reference's GSIS path is unusable as shipped.** `Init_Acceleration_New`
uses `AA_TMP` uninitialised; the six off-diagonal blocks it never writes hold
freed heap (values ~1e5) and corrupt every row of the HDG global matrix. The
shipped GSIS "converges" in 5 iterations to a field of magnitude 1e-27.
`tools/build_fortran.sh` zeroes it in its build copy. This port is correct by
construction.

**2. GSIS does not converge to the same discrete fixed point as CIS.** Both
schemes reach `residual_iterate < 1e-13`, but the *true* transport residual is
4.3e-14 for CIS and 3.2e-03 for GSIS: CIS is at the discrete kinetic fixed
point and GSIS is not. Their converged fields differ by 1.7e-2 at
`tau_R = 1e-1` — in the Fortran as well as here. The cause is structural: the
macroscopic system is discretised by HDG while the kinetic one is upwind DG,
and the gap does **not** close under mesh or order refinement (only as
`tau_R → 0`). Proposal 1 §10 lists agreement to `rtol=1e-8` as "the
load-bearing property for the whole benchmark"; **it does not hold**, and
anything built on this solver needs to grade against the true residual and the
analytic limit instead. `tools/fixed_point_study.py` reproduces the
measurement.

## Configuration

`cases/*.yaml` replaces the Fortran namelist one-to-one, plus the settings
that were implicit in the reference:

```yaml
iteration:  {tol: 1.0e-8, tmax: 8000000, true_residual: false}
scheme:     {accflag: 0, acc_variant: A, stabilisation: [1.0, 1.0, 1.0],
             scale_by_hmin: false, on_cycle: raise}
velmesh:    {npole: 20, nazim: 40}          # nazim is forced even
dg:         {deg: 3}
flow:       {cv: 1.0, vg: 1.0, tau_r: 1.0e-3, tau_n: 1.0e5, tau_thr: 1.0}
mesh:       {file: ../meshes/A1_Nx11_Ny11.msh}
boundaries:
  - {name: SWall, phyid: 11, type: thermalising, temp: 0.0}
  - {name: NWall, phyid: 12, type: thermalising, temp: 1.0}
  - {name: EWall, phyid: 13, type: thermalising, temp: 0.0}
  - {name: WWall, phyid: 14, type: thermalising, temp: 0.0}
restart:    {enabled: false, path: null, error_on_stale: true}
output:     {dir: ../out/run, field: true, runtime_log: true, run_record: true}
performance: {precompute_inverse: true, kernel: numba, threads: 1}
```

Boundary types are `thermalising`, `nonthermalising` (diffusely reflecting,
adiabatic) and `periodic` — all three are implemented and dispatched on, which
the reference does not do. `scheme.acc_variant: B` selects the `TAU_R`-rescaled
acceleration from `Synthetic_Acceleration1.f90`; it is ported faithfully, which
means it diverges (see the issues doc).

Every default reproduces the Fortran.

## Layout

```
pybte/
  quadrature.py velocity.py basis.py     angular + reference-element machinery
  mesh/                                  gmsh 2.2 reader, topology, geometry
  integration.py                         the precomputed integral tensors
  ordering.py                            per-direction topological sweep order
  sweep.py  _kernels.py                  the hot path (numba)
  sweep_reference.py                     the same sweep in plain numpy/LAPACK
  moments.py bc.py analytic.py io_output.py driver.py
  acceleration/                          GSIS; not imported at all when accflag=0
cases/  meshes/  docs/  tests/{unit,stage,integration}
tools/
  build_fortran.sh  run_fortran.py       build and drive the reference
  dump_fortran_stages.py compare_stages.py   the port's debugging instrument
  make_meshes.py                         mesh family (reproduces the shipped one)
  fixed_point_study.py validation_sweep.py
```

CIS and GSIS are separable at the module boundary: with `accflag: 0` nothing
in `pybte/acceleration/` is imported.

## Reproducing the Fortran side

```bash
tools/build_fortran.sh --variant A        # gfortran + MKL PARDISO
tools/build_fortran.sh --variant A --no-mkl   # reference LAPACK + dense stub
python tools/run_fortran.py --out /tmp/f90 --accflag 1 --dump
python tools/compare_stages.py /tmp/f90/dump
```

The reference tree is never modified. `build_fortran.sh` copies it, overlays
additive shims (MKL/OpenMP no-ops, an optional dense PARDISO stand-in, the
dump instrumentation) and applies the one correctness fix described above;
`--no-fixes` reproduces the broken build.

## Tests

```bash
pytest tests/unit                    # ~30 s, no Fortran needed
pytest tests/stage                   # needs fortran-reference/golden dumps
pytest -m "not slow"                 # everything quick
pytest                               # full, several minutes
```

`docs/`: `EQUATIONS.md` (model and discretisation), `INDEXING.md` (the 1-based
↔ 0-based mapping, in one place), `FORTRAN_ISSUES.md` (eleven findings).
