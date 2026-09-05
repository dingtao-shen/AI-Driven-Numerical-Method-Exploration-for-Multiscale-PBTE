# Proposal 2 — A benchmark for AI agents on stiff kinetic transport solvers

**Working title:** `StiffKinetic-Bench` (v0 instantiation: gray linear Callaway phonon BTE)
**Status:** plan of record
**Depends on:** Proposal 1 (`solver-python/` complete through Gate 6)
**Audience:** Claude Code, executing against this repository

---

## 1. Motivation

### 1.1 What exists

Benchmarks that ask an AI agent to write a PDE solver have appeared rapidly. PDEAgent-Bench
provides 645 instances across 6 mathematical categories and 11 PDE families with
library-specific tracks for DOLFINx, Firedrake and deal.II, each instance supplying an
agent-facing specification, a reference solution on a prescribed grid, and case-specific
accuracy and runtime targets, graded in stages: executability, then numerical accuracy,
then computational efficiency. PDE Agent Bench covers eight families — heat,
convection–diffusion, Stokes, Navier–Stokes, Helmholtz, biharmonic, linear elasticity,
reaction–diffusion — in 191 cases, each with an accuracy and time budget, executed in an
isolated environment against a per-case reference. Alongside these sit PDE-Bench, FEM-Bench,
CodePDE, AutoPDE and SciML-agents.

Two observations follow.

**First, they are all elliptic/parabolic/hyperbolic PDEs in the classical sense.** None
covers kinetic transport: no linear Boltzmann, no discrete ordinates, no neutron transport,
no radiative transfer, no rarefied gas, no phonon BTE. The entire class of equations whose
numerical difficulty lives in a stiff collision operator is absent.

**Second, and more importantly, they all grade the same two things: is the answer accurate,
and is it fast.** That is the correct grading for classical PDEs, where a competent
discretisation plus a black-box linear solver gets you there. It is the wrong grading for
kinetic transport, where the discretisation is often the easy part and the entire
difficulty is *iterative*.

### 1.2 The gap we can fill

For the linear Boltzmann family, source iteration is a Richardson iteration that isolates
the transport operator so it can be inverted by a sweep. It works well in many settings but
deteriorates in diffusive regimes; in optically thick, highly scattering settings the
contraction factor approaches unity, producing slow convergence or pseudo-convergence.
Acceleration is therefore not an optimisation, it is a precondition for the method being
usable at all. And getting acceleration right is delicate: on polytopic discretisations,
MIP-based DSA remains robust across optical thickness, scattering ratio, angular
quadrature, mesh refinement, polynomial degree and mesh anisotropy, while SIP-based DSA can
lose robustness in the intermediate regime.

The same story holds for phonon transport. Synthetic iterative schemes that tightly couple
macroscopic moment equations to the kinetic equation deliver one to three orders of
magnitude faster convergence than conventional implicit DOM in the near-diffusive regime,
and on non-gray Callaway benchmarks give speed-up factors in the several tens for both
iteration count and wall-clock.

This gives us a grading signal that no existing agent benchmark uses:

- **Hardware-independent.** Iteration counts, sweep counts, and measured contraction
  factors are identical on any machine. Existing benchmarks gate on wall-clock, which is
  fragile in containers and forces conservative thresholds. We can be strict.
- **Continuously parameterised by a physical stiffness knob.** Sweeping the Knudsen number
  over four decades gives graded difficulty and partial credit for free, rather than a
  single pass/fail.
- **Very hard to fake.** An agent cannot achieve bounded iteration counts across four
  decades of stiffness without actually implementing a correct acceleration scheme. And we
  can add a fixed-point-preservation check — the accelerated solution must equal the
  unaccelerated converged solution — which rules out the dominant cheat of loosening the
  tolerance or perturbing the physics.
- **Backed by theory.** Contraction-factor bounds and asymptotic-preservation are
  properties with known targets, not vibes.

### 1.3 One-line pitch

> Existing benchmarks ask whether an agent can *solve the equation*. We ask whether it can
> *build a solver that stays robust in the stiff limit* — the thing that actually separates
> a working transport code from a toy one.

### 1.4 Why this project, from this code

We hold an unusual asset: a research solver in which conventional source iteration (CIS)
and a general synthetic iterative scheme (GSIS) are the same codebase behind a single flag,
converging to the same discrete fixed point. That is precisely the `solver/` and `oracle/`
pair a benchmark task needs, and it exists because it was built for real research, not for
a benchmark. The parameter that controls difficulty (`TAU_R`) is already a config field.
The metric (iteration count) is already logged. The diffusion-limit reference (an analytic
Fourier series for the shipped cavity) is already implemented.

### 1.5 Scope honesty

Gray linear Callaway phonon BTE alone is too narrow for a standalone benchmark paper, and
reviewers will say so. The project is therefore framed and built as **stiff kinetic
transport**, with the phonon BTE as the first instantiation and a plug-in interface (§5.6)
so a second kinetic equation — neutron transport S_N, or a BGK model — can be added by a
collaborator without restructuring. The v0 paper may ship phonon-only, but the abstraction
must be in place from the start or it will never be retrofitted.

---

## 2. The idea, concretely

A task is a container plus an automatic grader. Inside the container is a working but
deliberately incomplete solver. The agent gets a terminal and a task description. When it
stops, a script decides pass or fail.

The v0 flagship task, in full:

- **`environment/`** — the Python solver with the entire `acceleration/` subpackage removed
  and the GSIS branch stripped from the driver. What remains is a correct, competent,
  unaccelerated CIS solver, plus five case files sweeping `TAU_R ∈ {1e0, 1e-1, 1e-2, 1e-3,
  1e-4}`, plus stored reference solutions for all five.
- **`task.md`** — "The solver converges for the first two cases and stalls on the rest.
  Modify it so all five converge to `tol=1e-8` in under 200 outer iterations each, without
  changing the converged solution."
- **`verifier/`** — runs all five cases and checks four things:
  1. all five report converged
  2. iteration count at `TAU_R=1e-4` is below 200
  3. converged temperature field matches the stored reference to `rtol=1e-6`
  4. the fixed point is unchanged — the solution equals the one the *unmodified* solver
     reaches on the two cases where it does converge, to `rtol=1e-8`
- **`oracle/`** — the GSIS implementation, restored. Must score 1.0, five times out of five.
- **null baseline** — the unmodified environment must score 0.0, five times out of five.

Check 4 is the load-bearing one. Without it an agent passes by raising the tolerance, by
damping the physics, or by declaring convergence on a different quantity. With it, the only
way through is to actually accelerate the iteration.

---

## 3. Solver assessment — is the Fortran code sufficient?

**Yes, and unusually cleanly.** Recorded here because it is the premise of everything below.

### 3.1 The CIS/GSIS separation is already at module boundaries

The driver loop is:

```fortran
DO
   CALL DG_Solver_VDF ()                              ! shared
   IF (ACCFLAG.EQ.0) CALL Calculate_Macro_Properties ()
   IF (ACCFLAG.EQ.1) THEN
      CALL Calculate_SRC_ACC_HoTfromDVM ()
      CALL Global_Problem_Solver_ACC ()
      CALL Local_Problem_Solver_ACC ()
      CALL Correct_VDF_Calculate_Macro_Properties ()
   END IF
   CALL Calculate_Residual_T (RESIDUL)                ! shared
   IF (RESIDUL.LT.TOL) EXIT
END DO
```

The acceleration lives entirely in one module, touched from exactly six call sites. Removing
it yields a complete, correct CIS solver with no dangling references. This is the ideal
starting shape for task construction and is rarer than it sounds.

### 3.2 Everything the benchmark needs is present

| Requirement | Present? | Where |
|---|---|---|
| Baseline that is correct but slow | Yes | `ACCFLAG=0` |
| Oracle that fixes it | Yes | `ACCFLAG=1`, `Synthetic_Acceleration.f90` |
| Both converge to same fixed point | By construction; verified at Proposal 1 Gate 4 | — |
| Stiffness knob | Yes | `TAU_R` in `control.in`; `Kn = TAU_R` for `Vg=L=1` |
| Second, independent stiffness knob | Yes | `TAU_N` — hydrodynamic limit |
| Iteration count as primary metric | Yes, already logged | `RunTime.txt` |
| Diffusion-limit analytic reference | Yes | 200-term Fourier series in `Out_Put_Result.f90` |
| Mesh/order/angle refinement axes | Yes | `DEG`, `NPOLE`, `NAZIM`, `.msh` |
| Realistic latent bugs for debug tasks | Yes — see §3.3 | — |

### 3.3 Gaps, and why each is an asset rather than a problem

Four issues found in the Fortran. Each is a real defect *and* a ready-made task.

1. **The BC dispatch is commented out.** In `Solvers.f90` all three `IF (BC_TYP...)`
   branches are commented; the thermalising branch is applied unconditionally to every
   boundary face. The non-thermalising machinery (`Calculate_FLUX_WALL`) exists but is never
   wired in; periodic is dead despite `Spatial_Mesh.f90` doing the face pairing.
   → Task T5.

2. **Two acceleration variants differing by a `TAU_R` rescaling of the macroscopic system.**
   `Synthetic_Acceleration.f90` has `OO/TAU_R` in the momentum block;
   `Synthetic_Acceleration1.f90` multiplies through by `TAU_R` instead. Algebraically
   identical, numerically different as `TAU_R → 0`. → Task T7.

3. **The residual is a change-between-iterates, not a true residual.** With a contraction
   factor near 1, successive iterates are close while both are far from the fixed point.
   The solver will report convergence at `tol=1e-8` while being materially wrong. → Task T3,
   probably the most interesting task in the suite.

4. **`A_SOL` is rebuilt and LU-factorised every element, every direction, every iteration**,
   despite depending only on iteration-invariant data. → Task T8.

Two further items are pure engineering and are fixed during the port (Proposal 1 §7.3,
§7.4): the silent restart-file read, and the quadratic global assembly.

### 3.4 What is genuinely missing and must be built

- Meshes at more than one refinement level (only `A1_Nx11_Ny11.msh` ships)
- Any test suite
- Any MMS infrastructure
- Reference solutions in a machine-readable format
- Non-gray/spectral physics — **out of scope, and fine**; it would multiply the work without
  changing the scientific claim

**Verdict: proceed.** The code supplies both halves of at least ten task families. Its
defects are, unusually, the most valuable part of it.

---

## 4. Repository layout

```
stiffkinetic-bench/
├── PROPOSAL_1_PYTHON_PORT.md
├── PROPOSAL_2_BENCHMARK_PROJECT.md
├── fortran-reference/           # read-only original + golden logs + stage dumps
├── solver-python/               # Proposal 1 deliverable; the single source of truth
├── tasks/
│   ├── t01-diffusive-acceleration/
│   │   ├── task.toml            # metadata, budgets, tags
│   │   ├── task.md              # agent-facing prompt
│   │   ├── ablation.yaml        # declarative spec of what to strip (see §5.2)
│   │   ├── environment/         # GENERATED — do not hand-edit
│   │   ├── verifier/
│   │   │   ├── test.sh
│   │   │   └── checks.py
│   │   ├── oracle/
│   │   │   └── patch.diff       # restores the removed capability
│   │   ├── reference/           # *.npz reference solutions
│   │   └── Dockerfile
│   └── ...
├── tools/
│   ├── make_env.py              # solver-python + ablation.yaml -> environment/
│   ├── make_reference.py        # generate + cross-check reference data
│   ├── validate_task.py         # oracle 5/5 pass AND null 5/5 fail
│   ├── run_experiments.py       # the model x scaffold x k grid
│   └── analyse_transcripts.py   # failure-mode coding support
├── experiments/
│   ├── configs/
│   └── results/
├── docs/
└── paper/
```

### 4.1 The one rule that keeps 30 tasks maintainable

**`environment/` is generated, never hand-edited.** Each task declares an `ablation.yaml`
saying which modules, functions, and config paths to remove or alter; `tools/make_env.py`
applies it to the current `solver-python/`. When a bug is fixed in the solver, regenerate
all thirty environments with one command. Hand-edited environments will drift within weeks
and the drift will not be noticed until results are already collected.

`tools/make_env.py` must also, unconditionally:

- strip `.git`
- strip all comments and docstrings mentioning the removed capability
- remove unused imports and empty stub functions left by the ablation (a `pass`-bodied
  `macro_solve()` is a giant hint)
- purge the removed capability's name from README, config schemas, and CLI help
- insert the benchmark canary string
- emit a manifest listing every file and its hash

---

## 5. Porting the solver into the benchmark

### 5.1 The solver is a dependency, not a copy

`solver-python/` is developed and tested independently and versioned (`pybte==0.3.1`). Task
environments contain a *derived* copy at a pinned version, recorded in `task.toml`. This
keeps solver development and benchmark development from blocking each other.

### 5.2 `ablation.yaml`

```yaml
base_version: pybte==0.3.1
remove:
  packages:   [pybte/acceleration]
  functions:
    - pybte/driver.py::_gsis_step
  config_keys:
    - scheme.accflag
    - scheme.acc_variant
    - scheme.stabilisation
    - flow.tau_thr
rewrite:
  - file: pybte/driver.py
    strategy: drop_branch
    branch: "case.scheme.accflag == 1"
scrub_terms: [GSIS, synthetic, acceleration, HDG, trace, PARDISO, splu, UQ, U_TRACE]
keep_tests:  [tests/unit, tests/integration/test_cis_*.py]
drop_tests:  [tests/integration/test_gsis_*.py, tests/stage/test_acc_*.py]
cases:       [tau_r_1e0, tau_r_1e-1, tau_r_1e-2, tau_r_1e-3, tau_r_1e-4]
```

`scrub_terms` is a hard requirement, applied to file contents *and* filenames. An agent that
greps `acceleration` and finds a leftover comment has been handed the answer.

### 5.3 Container

Base `python:3.11-slim`. Preinstall every dependency at build time — numpy, scipy, numba,
meshio, pyyaml, pytest. Nothing installed at runtime; the container runs with networking
disabled during grading. No MKL, no PARDISO, no compiler toolchain beyond what numba needs.
Pin `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `PYTHONHASHSEED=0`.
Target image size under 1.5 GB.

### 5.4 Case-size discipline

**This is the constraint most likely to blow up the project.** Research cases take minutes
to hours; a benchmark case must fit a 10–15 minute *total* verification budget, inside which
an agent will have run the solver many times while iterating.

Shrink by: coarser mesh (target `N_TRIS` 100–300), fewer angles (`NPOLE=10, NAZIM=20`),
lower `DEG` (2 rather than 3), and a `TMAX` cap tuned so the stalling cases truncate rather
than run forever.

**The shrunk case must preserve the pathology.** After shrinking, verify explicitly that
CIS iteration count still grows by orders of magnitude as `TAU_R` decreases, and that GSIS
still stays flat. A shrink that accidentally makes CIS converge quickly destroys the task.
This check goes in `tools/validate_task.py` and must run for every task.

**Do this on the very first task, before building the other twenty-nine.**

### 5.5 Reference data

Generated by the *full* solver at the task's exact discretisation, stored as `.npz`
(`temp`, `qx`, `qy`, `temp_dofs`, `iterations`, `config_hash`), kept under a few MB per task.

Cross-checks before a reference is accepted:
- CIS and GSIS agree to `rtol=1e-8` (where CIS converges at all)
- in the diffusive limit, agrees with the analytic Fourier series to discretisation error
- mesh- and angle-refined runs show monotone convergence toward it
- ballistic limit sanity: `TAU_R → ∞` reproduces the known ballistic behaviour

Never ship a reference that has only been produced by one code path.

### 5.6 The plug-in interface for a second physics

Define, and make the phonon BTE implement, an abstract kinetic problem:

```python
class KineticProblem(Protocol):
    def directions(self) -> Directions: ...
    def sweep(self, f, source) -> np.ndarray: ...
    def moments(self, f) -> Moments: ...
    def scattering_source(self, moments) -> np.ndarray: ...
    def stiffness_parameter(self) -> float: ...
```

Verifiers and task scaffolding are written against this interface, not against `pybte`
internals. A collaborator adding neutron transport S_N implements the protocol and inherits
the whole task-construction toolchain. Build this in Phase B, when there is exactly one
implementation — retrofitting an abstraction over thirty finished tasks does not happen.

---

## 6. Implementation phases and acceptance criteria

### Phase A — one task, end to end (3 weeks)

Build **only** T01. Then:

1. `tools/validate_task.py` reports **oracle 5/5 pass**
2. same tool reports **null baseline 5/5 fail**
3. shrunk case verified to preserve the CIS/GSIS iteration-count gap (§5.4)
4. total verification wall-clock under 15 minutes
5. **run Claude Opus 5 + Claude Code, k=5, and record how many pass**

Criterion 5 decides the project.

- **5/5 pass** → the task is too easy. Do not proceed. Escalate difficulty (deeper `TAU_R`,
  add the `TAU_N` hydrodynamic axis, tighten the fixed-point check) and repeat Phase A.
- **0/5** → good, but check it is failing for the *right* reason. Read all five transcripts.
  If they fail on environment friction (import errors, missing data, unclear prompt) rather
  than on the physics, that is a broken task, not a hard one.
- **1–2/5** → ideal. Proceed.

**Do not skip this.** One task costs three weeks and validates or kills the entire
direction. Thirty tasks cost six months and validate nothing until the end.

### Phase B — abstraction and tooling (2 weeks)

`KineticProblem` protocol; `make_env.py`, `make_reference.py`, `validate_task.py`
production-ready; task template; CI running `validate_task.py` on every task on every push.

Gate: T01 rebuilt from the template via the toolchain, still passing all Phase A criteria.

### Phase C — task suite (8 weeks)

Build T02–T10 (§7) with parameter variants. Each task must clear the Phase A criteria
individually. Track a live dashboard of task count, domain coverage, oracle status, and
baseline solve rate.

Gate: 12+ task templates, 30+ instances, all validated; measured frontier solve rate
between 5% and 40% (below 5% the benchmark cannot discriminate; above 40% it is too easy).

### Phase D — experiments (4 weeks)

§8. Gate: full result grid collected with `k=5`, failure-mode coding complete.

### Phase E — write-up (4 weeks)

Workshop-length paper first. Gate: submission, plus a public repo with a DOI.

**Total: ~21 weeks.** A workshop-scale cut (12 tasks, 3 models × 2 scaffolds, `k=3`, one
ablation) is reachable in ~12 weeks by trimming Phase C and D.

---

## 7. Task construction plan

Ten families, all realisable from this codebase. `solver/` and `oracle/` given for each.

### T01 — Diffusive-limit acceleration `[flagship]`
Baseline: CIS only. Ask: bounded iterations across `TAU_R ∈ {1e0 … 1e-4}`.
Oracle: GSIS. Checks: converged; iterations < 200 at stiffest; field matches reference;
fixed point preserved.

### T02 — Hydrodynamic-limit acceleration
Sweep `TAU_N → 0` at fixed large `TAU_R` — a *different* stiff limit, in which a
temperature-only diffusion acceleration does not help; the heat-flux equation with its
stress closure is required. Oracle: the full 7-field GSIS. This is the task that
distinguishes the suite from the neutron-transport literature, where the second relaxation
time has no analogue.

### T03 — Pseudo-convergence `[flagship]`
Ship the solver with its as-written iterate-difference residual. At `TAU_R=1e-4` it reports
converged at `tol=1e-8` while the field is materially wrong. Ask: diagnose why the reported
convergence is false and fix the criterion. Checks: field matches reference to `rtol=1e-6`;
the reported residual is now a true residual; no case is "fixed" by simply tightening `tol`
by six decades (cap the iteration budget so brute force cannot pass).
Tests understanding rather than coding, and is very hard to solve by retrieval.

### T04 — Asymptotic preservation on an under-resolved mesh
Coarse mesh with `h >> Vg*TAU_R`, `TAU_R=1e-4`. Checks: L2 error against the ported analytic
Fourier series below tolerance, on the coarse mesh.

### T05 — Restore the diffusely reflecting wall
Genuinely unfinished code: `Calculate_FLUX_WALL` exists but is disconnected, and the BC
dispatch is commented out (§3.3-1). Checks: net normal heat flux through an adiabatic
boundary zero to `1e-10`; global energy balance closed; cross-plane case matches reference.
Maximally authentic — it is literally an open TODO from the real work.

### T06 — Stabilisation robustness
`ST(2)`, `ST(3)` fixed at 1.0 with `/Hmin` scaling commented out. On anisotropic and refined
meshes the contraction factor degrades. Ask: make it bounded across the mesh family.
Checks: contraction factor below threshold across all meshes × `DEG ∈ {1,2,3}` × `TAU_R`
values. Directly parallels the known MIP-versus-SIP DSA robustness result.

### T07 — Conditioning of the macroscopic system
Ship variant B (the `TAU_R`-multiplied form). At `TAU_R=1e-6` it loses accuracy or
conditioning. Ask: diagnose and fix. Checks: correct solution at `1e-6`; global matrix
condition number below threshold.

### T08 — Redundant factorisation `[performance]`
Ship the version that rebuilds and factorises `A_SOL` every iteration. Ask: reduce cost by
≥5× without changing iteration count or converged answer. Graded on `factorisation_count`
and `sweep_count` — hardware-independent — with wall-clock only as a loose timeout.

### T09 — Method of manufactured solutions
Ask: construct an MMS for the Callaway BTE and demonstrate design order. Checks: fitted
convergence slope within ±0.15 of `DEG+1` for `DEG ∈ {1,2,3}`; the manufactured source is
verified independently by symbolic differentiation in the verifier.

### T10 — Sweep-cycle robustness
Provide a mesh on which some directions induce a cyclic dependency, hanging the topological
ordering. Ask: detect and handle. Checks: terminates; correct solution; cycle-breaking
logged.

### Variant axes
Each family expands via `TAU_R` range, `TAU_N` range, `DEG`, mesh refinement,
`NPOLE/NAZIM`, geometry (cavity / cross-plane / porous), BC mix. Target 30–40 instances
from 10–12 templates. **Report template count and instance count separately** in the paper;
inflating one into the other is the most common dishonesty in this genre.

### Contamination control
- Grade on *properties*, never on reproducing a named method. "Bounded iterations under
  stiffness" admits many valid solutions; "implement GSIS" is one search away.
- Never name GSIS, synthetic acceleration, DSA, or HDG in `task.md`.
- Use geometry/BC/material combinations not appearing in any published figure.
- Canary string in every task, per the Terminal-Bench convention.
- Hold out reference data; ship only what the verifier needs.
- Track which tasks have close public analogues and report this honestly in a limitations
  section rather than hoping reviewers miss it.

---

## 8. Experiment plan

### 8.1 Main grid

Models (≥4, spanning ≥2 vendors plus ≥1 open-weights): e.g. Claude Opus 5, GPT-5.x,
Gemini 3 Pro, DeepSeek, Qwen.
Scaffolds (≥2): one minimal bash-loop scaffold with no domain knowledge (mini-swe-agent
class), one full-featured coding agent (OpenHands or Claude Code).
`k=5` per cell.

30 tasks × 4 models × 2 scaffolds × 5 = **1200 rollouts** for the main table.

Report per-cell pass rate with bootstrap CIs, plus per-family breakdown. Fixing the backbone
while varying the scaffold, and vice versa, is what separates model capability from
scaffolding — this is the standard protocol in this literature and reviewers will expect it.

### 8.2 Ablations, in priority order

**A1 — Domain hint (highest value).** Two versions of each task: bare, and with an explicit
hint pointing at moment-based acceleration. If pass rate jumps, the deficit is *knowledge*;
if it does not, the deficit is *implementation ability*. These are very different claims
about frontier models and this single number is the most quotable result in the paper.
Cost: 300 extra rollouts on a 10-task subset.

**A2 — pass@k, k=1…10.** Whether the suite is search-limited or capability-limited.
**This also decides whether the route-two RL-environment project is worth starting**: a
steep pass@k curve means a search process has something to find; a flat one means it does
not. Cost: 500 extra rollouts on a 10-task subset, best model only.

**A3 — Reasoning effort.** Low/medium/max on the best model. Reviewers will ask. Cost: 300.

**A4 — Stiffness scaling.** Pass rate as a function of `TAU_R` within T01/T02. Produces the
paper's most legible figure: capability degrading smoothly along a physical axis. Nearly
free — it comes out of the main grid.

**A5 — Solver-provided-tests.** With and without the CIS test suite in the environment.
Tests how much agents exploit available verification. Cost: 300.

Total including ablations: **~2600 rollouts.**

### 8.3 Failure-mode analysis — the core contribution

Hand-code 200 failed transcripts against a codebook developed from the first 30:

- implemented an inconsistent acceleration (fixed point moved)
- weakened the convergence criterion instead of the iteration
- damped or altered the physics to force convergence
- correct scheme, wrong boundary treatment for the macroscopic system
- did not recognise the diffusion limit at all
- recognised the need, could not implement the moment closure
- ran out of context / gave up
- environment friction (excluded from scientific conclusions, counted separately)

Two coders on a 50-transcript overlap; report Cohen's κ.

**This must be done by the domain expert.** An ML researcher reading these logs can only
report "the agent failed". Distinguishing "implemented a plausible but inconsistent
correction" from "implemented a consistent correction with a wrong boundary condition"
requires knowing what consistency means here. It is the entire reason this benchmark is
worth more than a generic one, and it is the section reviewers will remember. Budget 2–4
weeks and do not compress it.

### 8.4 Budget

API credits: rough order **$5k–15k**, dominated by long-horizon rollouts on frontier
models. Routes: researcher credit programmes at each vendor; affiliating with an existing
benchmark effort that already has sponsor credits; or cutting to the workshop scale (~600
rollouts) first.

Sandbox compute (Modal or Daytona) is a much smaller line item.

### 8.5 Reporting standards

- pass rate with bootstrap CIs, never a bare point estimate
- template count and instance count reported separately
- oracle pass rate and null-baseline fail rate for every task, in an appendix
- full task list with public-analogue contamination assessment
- all transcripts released
- Zenodo DOI; Apache 2.0

---

## 9. Risks

| Risk | Sign | Response |
|---|---|---|
| Tasks too easy | Phase A gives 5/5 | Escalate before scaling; add `TAU_N` axis, deepen stiffness, tighten fixed-point check |
| Tasks too hard | All models 0% everywhere | Add a graded difficulty ladder per family so partial credit discriminates |
| Shrinking destroys the pathology | CIS converges fast on the small case | §5.4 check in `validate_task.py`, run on every task, every CI push |
| Single-physics narrowness | Reviewer: "why should I care about phonons" | §5.6 abstraction from Phase B; recruit a neutron-transport or rarefied-gas collaborator during Phase C |
| Contamination | Agent names GSIS unprompted | Property-based grading; scrub terms; report honestly |
| Environment friction masquerading as difficulty | Failures cluster on imports/paths | Track friction as a separate failure category; fix and re-run |
| Benchmark saturates in a year | Next model release solves it | Design the stiffness ladder to extend by two more decades without new tasks |
| Effort sink | Six months, no paper | Workshop-first; ship at 12 tasks |

---

## 10. Immediate next actions

1. Execute Proposal 1 through Gate 4. Nothing here starts before CIS and GSIS are verified
   to reach the same fixed point.
2. Build **T01 only**, with the shrunk case, and run `validate_task.py`.
3. Run Claude Opus 5 + Claude Code, `k=5`, on T01. Record the number.
4. Read all five transcripts, regardless of the score.
5. Decide, on that evidence, whether to proceed to Phase B, escalate difficulty and repeat
   Phase A, or stop.

Step 5 is a real decision point, not a formality.
