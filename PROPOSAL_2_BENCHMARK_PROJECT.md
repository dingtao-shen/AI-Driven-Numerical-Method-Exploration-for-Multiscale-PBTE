# Proposal 2 — A benchmark for AI agents on stiff kinetic transport solvers

**Working title:** `StiffKinetic-Bench` (v0 instantiation: gray linear Callaway phonon BTE)
**Status:** plan of record — v0.2 (2026-09-08): grading is gates + score; the
"requires domain knowledge" claim of v0.1 was tested and withdrawn (§1.2, §6)
**Depends on:** `solver-python/` — complete, validated, 156 tests green
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
- **Hard to fake — not hard to solve.** A fixed-point-preservation check (the accelerated
  solution must equal the unaccelerated converged solution, certified by a transport
  residual computed outside the submission) rules out the dominant cheats: loosening the
  tolerance, damping the physics, declaring convergence on another quantity. It does *not*
  make the task require domain knowledge. Full GMRES on the source-iteration operator, one
  sweep per product, passes every shipped cell, and a frontier agent found it unaided in
  44 minutes (§6). The benchmark's signal is therefore the **score spread** across agents,
  models and scaffolds and the gate pass rate of weaker ones — not whether the strongest
  model passes. v0.1 claimed otherwise; the claim was measured and withdrawn.
- **Backed by theory.** Contraction-factor bounds and asymptotic-preservation are
  properties with known targets, not vibes.

### 1.3 One-line pitch

> Existing benchmarks ask whether an agent can *solve the equation*. We ask whether it can
> *build a solver that stays robust in the stiff limit* — the thing that actually separates
> a working transport code from a toy one.

### 1.4 Why this project, from this code

We hold a research solver in which source iteration (CIS) and a published synthetic
acceleration (GSIS) are one codebase behind a flag, validated bit-for-bit against its
Fortran original and against the publication's own 2-D table (unaccelerated counts 5/5
identical, accelerated 4/5). The stiffness knobs (`tau_R`, `tau_N`) are config fields; the
metric (outer iterations = sweeps) is logged; the diffusion-limit series is implemented.

Two facts, both measured, shape what the benchmark can honestly claim:

**The published accelerated scheme does not converge to source iteration's fixed point.**
Its last step blends the synthetic solution with the kinetic moment under a damping factor,
and at a fixed point that blend leaves a stationary offset. On the publication's own five
test points, at its own discretisation, the displacement is `3.7e-3 .. 2.3e-2`; it is
present at every Knudsen number and every boundary type tested and does not vanish under
refinement (`solver-python/docs/LIMITATIONS.md` #1). A grader that requires agreement with
the unaccelerated fixed point therefore rejects the published scheme. That is a result,
and the grader is what produces it.

**A standard Krylov method does converge to it.** GMRES on `(I - T) u = g`, each sweep one
matrix–vector product, reaches source iteration's fixed point to round-off (transport
residual `1e-15`, agreement with CIS to `4e-11 .. 1.6e-7` wherever CIS can reach the
answer) in 228 sweeps where CIS needs 16 830 and in 992 where CIS does not converge in
200 000. It is the oracle (`scheme.method: krylov`, `solver-python/docs/KRYLOV.md`). It
is not Knudsen-independent, its cost is memory, and it does not apply to periodic faces
— which is where the next family (§7, F3) comes from.

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
deliberately unaccelerated solver. The agent gets a terminal and a task description. When
it stops, a script decides whether every gate passed and, if so, how well.

The v0 task, `tasks/t01-square`, in full:

- **`environment/`** — the Python solver with every acceleration path removed (the
  moment-based subpackage, the Krylov driver, their config keys), scrubbed of every term
  naming either, canaried. Plus ten case files: two boundary families on the same square
  and the same shrunk discretisation (200 elements, `DEG = 2`, `10 x 20` angles,
  `tol = 1e-8`), each at the publication's five `(Kn_R, Kn_N)` pairs.
- **`task.md`** — "Make every case converge — to the same answer the solver already
  converges to where it can — and make it fast." States the gates and the score rule.
- **`verifier/`** — owns its own copies of the cases, the mesh and a pristine solver.
  **Gates**, every cell: (1) converged within cap; (2) transport residual
  `||A f - b||/||b|| < 1e-7`, computed from the submitted distribution in the pristine tree;
  (3) field within `1e-5` of the certified reference, both as reported and as recomputed
  from the distribution. **Score**, only if every gate of every family passes: the
  speed-up in outer iterations over the unaccelerated solver on the cells where it is
  slow, geometric mean per family and then across families; censored where the
  unaccelerated solver never converged within its 200 000 calibration cap. Peak memory and
  wall-clock per iteration relative to the unaccelerated solver are reported alongside.
- **`oracle/`** — the Krylov driver restored and made the default. Must pass every gate.
- **null baseline** — the environment unmodified. Must fail.

| cell | regime | unaccelerated | oracle | speed-up |
|---|---|---|---|---|
| F1 (0.001, 1e5) | deep diffusive | > 200 000 | 992 | >= 201x |
| F1 (0.01, 1e5) | diffusive | 16 830 | 228 | 74x |
| F1 (0.1, 1e5) | transition | 317 | 37 | gate only |
| F1 (1, 1) | ballistic | 33 | 21 | gate only |
| F1 (10, 0.01) | hydrodynamic | 2 706 | 574 | 4.7x |
| F2 (0.001, 1e5) | deep diffusive | > 200 000 | 1 104 | >= 181x |
| F2 (0.01, 1e5) | diffusive | 32 343 | 281 | 115x |
| F2 (0.1, 1e5) | transition | 643 | 47 | gate only |
| F2 (1, 1) | ballistic | 75 | 23 | gate only |
| F2 (10, 0.01) | hydrodynamic | 2 444 | 688 | 3.6x |

F1: all walls isothermal, `T = 1` north. F2: east and west walls diffusely reflecting.
The oracle scores 42x. Validation: oracle gate PASS 5/5 (~40 s each), null FAIL 5/5
(~330 s), pathology preserved in both families, worst verification 346 s of a 900 s budget.

Gates 2 and 3 are the load-bearing pair. Gate 2 alone is not sufficient at small `Kn`
— the operator is ill-conditioned there and a `4e-12` residual was measured next to a
`1.5e-5` error — which is why gate 3 exists and why the references are cross-checked
against source iteration wherever it converges.

**What the score is for.** A method that passes every gate can still be slow, memory-hungry
or regime-dependent; the score, per family, is where that shows. The (10, 0.01) cells are
the Krylov oracle's weakest at 4–5x, and a moment-preconditioned method would be expected
to do better there; the score is what would reward it.

## 3. Solver assessment

**Sufficient, and unusually cleanly so.** Recorded here because it is the premise of
everything below.

### 3.1 The CIS/GSIS separation is already at module boundaries

The driver loop is:

```python
while True:
    self.sweep()                      # shared
    if self.acc is None:
        compute_moments(...)          # CIS
    else:
        self.acc.apply()              # GSIS: hot source, global solve, local solve, blend
    res = residual_iterate(...)       # shared
    if res < tol: break
```

The whole acceleration is one subpackage, `pybte/acceleration/`, imported lazily and touched
from a single call site. With `scheme.accflag: 0` it is never imported at all, so removing
it yields a complete, correct CIS solver with no dangling references. This is the ideal
starting shape for task construction and is rarer than it sounds.

### 3.2 Everything the benchmark needs is present

| Requirement | Present? | Where |
|---|---|---|
| Baseline that is correct but slow | Yes | `scheme.accflag: 0` |
| Oracle that fixes it | Yes | `scheme.accflag: 1` |
| Oracle that preserves the CIS fixed point | Yes — full GMRES on the outer iteration; **not** the shipped GSIS (§1.4) | `scheme.method: krylov`, `pybte/krylov.py` |
| Certificate of correctness | Transport residual (necessary) + certified reference (sufficient at small `Kn`) | verifier gates 2 and 3 |
| Stiffness knob | Yes | `flow.tau_r`; `Kn = TAU_R` for `Vg=L=1` |
| Second, independent stiffness knob | Yes | `TAU_N` — hydrodynamic limit |
| Iteration count as primary metric | Yes | `RunRecord.iterations`, `.sweep_count` |
| Diffusion-limit analytic reference | Yes | `pybte.analytic` |
| Mesh/order/angle refinement axes | Yes | `DEG`, `NPOLE`, `NAZIM`, `.msh` |
| Live pathologies for debug tasks | Yes — see §3.3 | — |

### 3.3 Task material already in the solver

Two kinds, and the distinction matters for how `ablation.yaml` is written.

**Live pathologies** — present in the working solver, no ablation needed. These are the
most valuable tasks because the difficulty is real rather than manufactured.

1. **The residual is a change between iterates, not a true residual.** With a contraction
   factor near 1 successive iterates are close while both are far from the fixed point.
   Measured at `TAU_R = 1e-2`, `tol = 1e-4`: 2 599 iterations, reported residual `9.99e-05`,
   relative error in `int T dA` `1.28e-01` against the same scheme converged to
   `1e-12` — **three orders of magnitude larger**. At `TAU_R = 1e-4` a truncated run is 89%
   wrong in `int T dA` while reporting a residual of `2.5e-06`. → **T03**.
2. **The published accelerated scheme converges to a different discrete fixed point than
   the unaccelerated one** (§1.4). Under property-based grading this is not a separate
   task — replacing the scheme with Krylov passes — so it is reported as a *finding* and
   folded into T01's gates (see T11).
3. **Acceleration variant B diverges.** Algebraically the same system as variant A,
   rescaled by `TAU_R`. It diverges geometrically from the first iteration at *every*
   Knudsen number tested, reaching NaN within ~30 iterations, and not because of
   conditioning — the two global matrices have comparable 1-norm condition estimates.
   → **T07**.

**Capabilities to ablate** — correct in the solver, removed by `ablation.yaml` to create
the task. Not defects; the environment generator takes them out.

4. Both acceleration schemes (`scheme.accflag`) → T01, T02.
5. The boundary-condition dispatch, including the diffusely reflecting wall and the
   periodic pairing → T05.
6. `scheme.scale_by_hmin`, the `/Hmin` scaling of the HDG stabilisation → T06.
7. `performance.storage: stored`, which LU-factorises the per-direction operators once at
   setup. Shipping `onthefly` restores the rebuild-and-refactorise-every-iteration
   behaviour → T08.
8. `scheme.on_cycle`, the sweep-graph cycle detection → T10.

**Verdict: proceed.** The solver supplies both halves of at least eleven task families.

### 3.4 What is genuinely missing and must be built

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
├── PROPOSAL_2_BENCHMARK_PROJECT.md
├── solver-python/               # the single source of truth
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

**A reference is certified by its transport residual, never by an iterate tolerance.**
This matters more than it sounds. The oracle's outer loop inherits the same round-off floor
the acceleration does, so at `TAU_R = 1e-3` it cannot meet an *iterate* tolerance of
`1e-10` even after 60 000 iterations — while its *transport* residual sits at `1.2e-9`,
i.e. the state is converged and the criterion is lying about it. That is T03's lesson
applied to our own tooling, and getting it backwards would silently ship references worse
than the runs they grade.

Cross-checks before a reference is accepted:
- transport residual `||A f - b||/||b|| < 1e-8` (the Krylov oracle delivers `1e-15`). Necessary,
  not sufficient on its own at small `Kn`, hence the next line. The shipped GSIS leaves it at `1.6e-6` to `5.6e-3` and must never generate a
  reference.
- CIS and the fixed-point-preserving oracle (§1.4) agree to `rtol=1e-4` where CIS converges
  at all — `TAU_R >= 1e-2` on the shrunk case. Do **not** substitute the shipped GSIS
  here: it disagrees with CIS by up to `1.7e-2` by design.
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

### Phase A — one task, end to end

Build **only** the v0 task. Then:

1. `tools/validate_task.py` reports **oracle gate PASS k/k** and **null baseline FAIL k/k**
2. the shrunk case preserves the pathology (§5.4)
3. total verification wall-clock under 15 minutes
4. **run Claude Opus 5 + Claude Code, k=5, and record gate pass rate and scores**

Criterion 4 decides the next step — but not by pass rate alone. Under gates + score:

- **every trial passes and the scores cluster** → the task does not discriminate at the
  top; difficulty must come from a new family (§7, F3), not from tightening budgets;
- **every trial passes and the scores spread** (≥ 3x between trials or models) → the
  benchmark has a signal at the top; proceed and add weaker models and scaffolds;
- **0/5 pass** → read all five transcripts. Failures on imports, paths or an unclear
  prompt mean a broken task, not a hard one.

**Record so far.** The v0.1 precursor (isothermal family only, budget-graded): 1/1 solved
by Claude Opus 5 + Claude Code in 44 minutes and $7, with a self-written full GMRES that
outperformed the then-oracle 3–8x. That result is what retired the v0.1 claim. The v0.2
task `t01-square`, k=5 (2026-09-08): **4 valid trials, 4/4 pass, scores 989x – 3 293x**
under the sweep-count rule then in force (spread 3.3x); one trial invalid (account quota).
Every valid trial built a physics-based preconditioner for GMRES — three a coarse-angle
transport operator, one a moment system — reaching near-Knudsen-independent counts of
6–37 sweeps. The score rule was then changed to sweep-equivalents (§2); the solutions were
lost to a reboot before they could be re-scored. **Second k=5 (2026-09-10), final rule** — unit =
the submission's own single-threaded sweep, capped by the verifier's: **5/5 pass, 548x – 974x**
(geometric mean 762x, spread 1.8x; oracle 34x). Every trial built a low-order preconditioner for
GMRES (coarse-angle x3, sparse low-order x1, diffusion-type x1); none touched the sweep kernel,
which the unit makes pointless; 3/5 hit the 2 h cap with a passing solution in place. An
intermediate wall-clock-unit rule was tried and withdrawn after one trial spent its budget
compiling sweep kernels. `experiments/results/t01-square*/NOTES.md`.

### Phase B — abstraction and tooling (2 weeks)

`KineticProblem` protocol; `make_env.py`, `make_reference.py`, `validate_task.py`
production-ready; task template; CI running `validate_task.py` on every task on every push.

Gate: T01 rebuilt from the template via the toolchain, still passing all Phase A criteria.

### Phase C — task suite (8 weeks)

Build T02–T11 (§7) with parameter variants. Each task must clear the Phase A criteria
individually. Track a live dashboard of task count, domain coverage, oracle status, and
baseline solve rate.

Gate: 12+ task templates, 30+ instances, all validated; across the tested models the
score spread is at least 3x on every scored family and the gate pass rate is not saturated
for at least two models. A suite where every agent gets the same score is not a benchmark.

### Phase D — experiments (4 weeks)

§8. Gate: full result grid collected with `k=5`, failure-mode coding complete.

### Phase E — write-up (4 weeks)

Workshop-length paper first. Gate: submission, plus a public repo with a DOI.

**Total: ~21 weeks.** A workshop-scale cut (12 tasks, 3 models × 2 scaffolds, `k=3`, one
ablation) is reachable in ~12 weeks by trimming Phase C and D.

---

## 7. Task construction plan

Eleven families, all realisable from this codebase. `solver/` and `oracle/` given for each.

### T01 — Square-domain acceleration, two wall families `[flagship, built]`
`tasks/t01-square`, §2. Baseline: source iteration only. Ask: converge every cell to the
shipped solver's own fixed point, fast. Oracle: Krylov (§1.4). Gates and score as in §2.
Families: F1 isothermal, F2 diffusely reflecting side walls. **F3 (periodic side walls)
is next**: it is the one configuration where the unaided GMRES solution declares itself
inapplicable (the partner element cannot be placed upwind, so one sweep is no longer an
exact evaluation of the map) and falls back to source iteration — measured: it fails the
convergence gate on the two stiff cells. F3 needs a periodic-capable oracle first; the
periodic temperature-jump boundary of the publication's long-film test is not implemented.

### T02 — Hydrodynamic-limit acceleration
Sweep `TAU_N → 0` at fixed large `TAU_R` — a *different* stiff limit, in which a
temperature-only diffusion acceleration does not help; the heat-flux equation with its
stress closure is required. Oracle: the full 7-field GSIS. This is the task that
distinguishes the suite from the neutron-transport literature, where the second relaxation
time has no analogue.

### T03 — Pseudo-convergence `[flagship]`
Ship the solver with its as-written iterate-difference residual. Measured: at
`TAU_R = 1e-2`, `tol = 1e-4` it stops after 2 599 iterations reporting `9.99e-05` while the
relative error in `int T dA` is `1.28e-01` (reference: the same scheme at `tol = 1e-12`)
— **three orders of magnitude larger**; at `TAU_R = 1e-4` a truncated run is 89%
wrong in `int T dA` while reporting `2.5e-06`. Ask: diagnose why the reported convergence
is false and fix the criterion. Checks: field matches reference to `rtol=1e-6`; the
reported residual is a true residual; no case is "fixed" by tightening `tol` alone (cap the
iteration budget so brute force cannot pass).
Tests understanding rather than coding, and is very hard to solve by retrieval.

### T04 — Asymptotic preservation on an under-resolved mesh
Coarse mesh with `h >> Vg*TAU_R`, `TAU_R=1e-4`. Checks: L2 error against the ported analytic
Fourier series below tolerance, on the coarse mesh.

### T05 — Restore the diffusely reflecting wall
Ablate the boundary-condition dispatch down to the thermalising branch, leaving the
diffusely reflecting wall and the periodic pairing unreachable. Checks: net normal heat flux
through an adiabatic boundary zero to `1e-10`; global energy balance closed; cross-plane
case matches reference. The tangential projection of the trace heat flux is the subtle part
— without it the accelerated run diverges rather than merely erring.

### T06 — Stabilisation robustness
Ablate `scheme.scale_by_hmin`, leaving the HDG stabilisation fixed at 1.0. On anisotropic
and refined meshes the contraction factor degrades. Ask: make it bounded across the mesh
family.
Checks: contraction factor below threshold across all meshes × `DEG ∈ {1,2,3}` × `TAU_R`
values. Directly parallels the known MIP-versus-SIP DSA robustness result.

### T07 — A divergent acceleration variant
Ship variant B, the `TAU_R`-multiplied form of the macroscopic system. It is algebraically
identical to variant A and diverges geometrically from the first iteration at **every**
Knudsen number tested, reaching NaN in ~30 iterations. It is not a conditioning problem —
the two global matrices have comparable 1-norm condition estimates — so the obvious
diagnosis is a false lead. Ask: diagnose and fix. Checks: converges at every `TAU_R` in the
ladder; transport residual at round-off; field matches reference.

### T08 — Redundant factorisation `[performance]`
Ship `performance.storage: onthefly`, which rebuilds and LU-factorises the per-direction
operators every element, every direction, every iteration, despite their depending only on
iteration-invariant data. Ask: reduce cost by ≥5× without changing iteration count or
converged answer. Graded on `factorisation_count`
and `sweep_count` — hardware-independent — with wall-clock only as a loose timeout.

### T09 — Method of manufactured solutions
Ask: construct an MMS for the Callaway BTE and demonstrate design order. Checks: fitted
convergence slope within ±0.15 of `DEG+1` for `DEG ∈ {1,2,3}`; the manufactured source is
verified independently by symbolic differentiation in the verifier.

### T10 — Sweep-cycle robustness
Provide a mesh on which some directions induce a cyclic dependency, hanging the topological
ordering. Ask: detect and handle. Checks: terminates; correct solution; cycle-breaking
logged.

### T11 — Fixed-point consistency of the acceleration `[withdrawn as a task]`
Shipping the published accelerated solver and asking for its fixed point to be repaired
without losing speed is, under property grading, solved by discarding the scheme and using
Krylov — the same solution as T01. §7's own rule forbids grading on the method, so T11 is
not a separate task. What it targeted survives as T01's gate 3 and as a reported finding:
the published scheme fails that gate on every published test point (§1.4).

### Variant axes
Each family expands via `TAU_R` range, `TAU_N` range, `DEG`, mesh refinement,
`NPOLE/NAZIM`, geometry (cavity / cross-plane / porous), BC mix. Target 30–40 instances
from 11–13 templates. **Report template count and instance count separately** in the paper;
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

1. Record the `t01-square` k=5 result (gate pass rate, score range) in §6 and in the root
   README; read every transcript, code the approaches.
2. **F3, periodic side walls**: implement the periodic temperature-jump boundary; extend the
   Krylov oracle to periodic faces (carry the periodic-face outflow in the state, or close
   each direction's periodic coupling exactly); calibrate, certify references, validate;
   add to `t01-square` as a third family and re-run k=5. This is the one family with
   evidence of discriminating against the unaided GMRES solution.
3. Tier-1 families that are not about the outer iteration and therefore not Krylov-solvable:
   T05 (boundary-condition restoration), T09 (manufactured solutions), T10 (sweep cycles).
4. Docker isolation (`run_rollouts.py --isolation docker`) before any number is published:
   directory isolation leaves the full solver on the same filesystem.
5. A second model and a second scaffold on `t01-square`, to measure the spread the score
   is designed to expose.
