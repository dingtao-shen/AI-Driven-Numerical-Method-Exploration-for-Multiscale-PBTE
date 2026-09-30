> Numerical backend of AI-Driven Numerical Method Exploration for Multiscale PBTE.
> Historical measurements are retained below; current scope and test evidence are in
> [project status](../docs/STATUS.md). Research contracts remain pending.

# pybte

Steady 2D2V gray linear **Callaway phonon Boltzmann** solver: nodal
discontinuous Galerkin in space, discrete ordinates in angle, with two
iteration paths —

* **CIS** — conventional source iteration (`scheme.accflag: 0`)
* **GSIS** — general synthetic iterative scheme, HDG macroscopic acceleration
  (`scheme.accflag: 1`)
* **Krylov** — GMRES on the source-iteration system, including reflecting and periodic state
  (`scheme.method: krylov`)

Ported from a Fortran research solver (`ACC_2D2V_LinearCallawayModel`) and
verified against it before that reference was retired; see `VALIDATION.md`.
Self-contained: no Intel compiler, no MKL, no PARDISO — `pip install -e .` and
PyPI wheels only.

## Install and run

```bash
pip install -e ".[dev]"          # numba + pytest; plain install works without
python -m pybte info cases/cavity_tauR1e-3_cis.yaml
python -m pybte run  cases/cavity_tauR1e-3_gsis.yaml
python -m pybte convert path/to/control.in     # old Fortran namelist -> yaml
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
record.residual_history    # (n_iter,)  the original iterate residual
record.residual_true       # (n_iter,) or None -- a true transport residual
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

Historical measurements below concern the shipped implementation and sampled cases, not a
new-domain robustness claim. Current engineering test results are in ../docs/STATUS.md.

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

### Do not trust a CIS residual

The stopping criterion measures the *step* between iterates, not the *error*.
Source iteration converges linearly, so the error left when you stop is
`r * rho/(1-rho)`, and `rho -> 1` as the medium becomes optically thick:

| `tau_R` | CIS iterations | reported residual | actual error in `int T dA` |
|---|---|---|---|
| 1e-2 | 16 836 (converged) | 1.0e-08 | 1.0e-05 |
| 1e-3 | 200 000 (truncated) | 1.5e-06 | **4.4e-02** |
| 1e-4 | 200 000 (truncated) | 2.5e-06 | **2.2e-01** — 89% wrong |

`record.error_estimate` fits `rho` to the tail of the residual history and
reports the implied error; it lands within an order of magnitude of the truth.
`record.contraction` gives `rho` itself.

## Fidelity

Verified against the Fortran reference stage by stage, on every iteration,
before that reference was removed from the tree (it remains in git history at
the `Port ACC_2D2V_LinearCallawayModel` commit):

* **Setup was bit-identical.** Mesh topology, geometry, angular quadrature,
  both nodal bases, all nine integral tensors and the sweep ordering
  reproduced the Fortran exactly — not merely to tolerance. That required
  matching the reference's summation orders and its use of libm `pow` for
  monomials.
* **CIS was bit-identical over full runs**, the longest verified being 16 836
  iterations at `tau_R = 1e-2`: identical iteration count, and every entry of
  the residual and mass histories equal to the last bit. The element LU is a
  hand-written `DGETF2`/`DGETRS` look-alike for exactly this reason, which
  also let the factorisation be hoisted out of the iteration.
* **GSIS agreed to ~2.5e-13** relative, stage by stage; the assembled global
  matrix matched the Fortran CSR to 1.5e-15 with identical sparsity. The
  residual is the sparse LU reordering the trace solve relative to PARDISO.
* **The sampled output field matched** the reference's own Tecplot file to
  `5e-8` — exact at that file's `ES16.6` print precision.

`VALIDATION.md` has the full sweep. Those results are frozen as assertions in
`tests/integration/test_regression.py`, so the solver cannot drift away from
them now that the reference is gone.

## Two things you should know before trusting a result

Both are documented in full in **`docs/LIMITATIONS.md`**, along with two more.
The accelerated path that keeps CIS's fixed point is `scheme.method: krylov`
(`docs/KRYLOV.md`); `docs/DEFECT_CORRECTION.md` documents an earlier experimental correction.

**1. The current implementation has recorded CIS/GSIS fixed-point differences.**
The retained measurements include true transport residuals 4.3e-14 (CIS) and
3.2e-03 (GSIS), and a field gap 1.7e-2 at `tau_R = 1e-1`.
These observations do not establish a universal property of GSIS or rule out
an implementation issue. The causal interpretation and discrete compatibility
need further verification. `tools/fixed_point_study.py` and its JSON data are
preserved; no numerical formula was changed in phase0.
`scheme.defect_omega` remains an experimental option (`docs/DEFECT_CORRECTION.md`).

**2. Do not trust a CIS residual.** The stopping criterion measures the *step*
between iterates, not the *error*, and source iteration's contraction factor
approaches 1 as the medium becomes optically thick. At `tau_R = 1e-4` CIS
reports a residual of `2.5e-6` after 200 000 iterations while its answer is
89% wrong. `record.error_estimate` reports the implied error; look at it
before believing a CIS result. Table below.

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
the reference does not do. `scheme.acc_variant: B` selects the
`TAU_R`-rescaled acceleration; it is ported faithfully, which means it
diverges (`LIMITATIONS.md` #3).

Every default reproduces the original solver.

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
cases/  meshes/  docs/  tests/{unit,integration}
tools/
  make_meshes.py        mesh family (reproduces the original shipped mesh)
  benchmark.py          the performance targets
  validation_sweep.py   regenerates VALIDATION.md
  fixed_point_study.py  measures the CIS/GSIS gap and how it scales
```

CIS and GSIS are separable at the module boundary: with `accflag: 0` nothing
in `pybte/acceleration/` is imported.

## Tests

```bash
pytest -m "not slow"                 # ~60 s
pytest                               # full, about 20 minutes
```

The slow half is `tests/integration/test_regression.py`, which re-runs the
converged cases and checks their iteration counts and masses against frozen
values. Several of those are CIS runs of 16 000+ iterations.

`docs/`: **`LIMITATIONS.md`** (read this one), `KRYLOV.md` (GMRES on the
outer iteration: same fixed point as CIS, `scheme.method: krylov`),
`DEFECT_CORRECTION.md` (experimental repair of the GSIS fixed point, off by
default),
`EQUATIONS.md` (model and
discretisation), `INDEXING.md` (the array-layout conventions).
