# Archive

A benchmark for AI coding agents on a problem class no existing PDE-agent
benchmark covers: **iterative robustness in stiff kinetic transport**.

The agent is given a correct but unaccelerated solver for the steady phonon
Boltzmann equation with dual relaxation times (Callaway model) — nodal
discontinuous Galerkin in space, discrete ordinates in angle, source
iteration in time. Source iteration converges in a handful of steps where
the medium is thin and not at all where it is optically thick or
hydrodynamic. The task is to deliver an accelerated solver that reaches
**the same discrete fixed point** on every case, and it is scored on how much
faster it gets there.

## What is claimed, and what is not

* **Grading is on properties of the delivered solver, never on the method.**
  Correctness is decided by a transport residual computed from the
  submitted distribution in a pristine tree, and by agreement with a
  certified reference; speed is the speed-up over source iteration in
  sweep-equivalents, timed by the grader. Loosening a tolerance, damping the physics or declaring
  convergence on another quantity cannot pass.
* **Frontier agents pass, with domain knowledge.** Full GMRES on the
  source-iteration operator already passes every gate (oracle: 27x on the
  three-family task); in three k=5 runs every valid Claude Opus 5 trial
  went further and built a physics-based low-order operator — as a GMRES
  preconditioner, or once as a classical synthetic correction — that makes
  the iteration count nearly Knudsen-independent, in 1–2 hours each,
  scoring 456x – 1 157x on the three-family task under the final rule. The
  benchmark's signal is therefore the *score spread* across agents, models
  and scaffolds (2.5x within one model here, and it separates the
  preconditioners), and the gate pass rate of weaker ones — not whether the
  strongest model passes.
* **The published accelerated scheme does not pass.** The solver reproduces
  the published GSIS iteration counts on the published 2-D test points
  exactly (CIS 5/5, GSIS 4/5), and at every one of those points GSIS
  converges to a fixed point displaced from source iteration's by
  `3.7e-3 .. 2.3e-2` — a property of its damped blend, not a bug. The grader
  detects it. `solver-python/docs/LIMITATIONS.md` #1.

## Layout

```
solver-python/          the solver (pybte): CIS, GSIS, Krylov; 156 tests
tasks/t01-square/       the task: three wall families x five (Kn_R, Kn_N) cells
  ablation.yaml         what make_env.py removes to build environment/
  environment/          GENERATED - what the agent sees; scrubbed, canaried
  verifier/             gates + score; its own pristine solver, cases, mesh
  reference/            certified reference solutions (*.npz)
  oracle/               a passing solution, as a patch against environment/
tools/                  make_env, calibrate_cells, make_reference, make_oracle,
                        validate_task, run_rollouts
experiments/results/    every rollout: transcript, per-cell scores, summary
PROPOSAL_2_BENCHMARK_PROJECT.md   design and plan
```

## The task

Same square, same discretisation (200 elements, `DEG = 2`, `10 x 20`
angles, `tol = 1e-8`), three boundary families, the five `(Kn_R, Kn_N)` pairs
of the reference publication's 2-D test:

| cell | regime | source iteration | oracle (Krylov) | speed-up (sweep count) |
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
| F3 (0.001, 1e5) | deep diffusive | > 200 000 | 1 063 | >= 188x |
| F3 (0.01, 1e5) | diffusive | 32 230 | 286 | 113x |
| F3 (0.1, 1e5) | transition | 624 | 55 | gate only |
| F3 (1, 1) | ballistic | 70 | 27 | gate only |
| F3 (10, 0.01) | hydrodynamic | 1 920 | 775 | 2.5x |

F1: all four walls isothermal. F2: east and west walls diffusely
reflecting (adiabatic). F3: east and west walls periodic — the one family
where a sweep is no longer an exact function of the moments alone, because
a periodic face's partner cannot always be placed upwind; the oracle carries
the partner's outflow in its state. Score = geometric mean of the speed-ups per family,
then across families, in **sweep-equivalents** (end-to-end wall-clock over
the time of one fine sweep, both timed by the verifier in-session, setup
included) so that work moved out of the sweep is paid for at its real cost.
The oracle scores 27x under this rule (42x by raw sweep count, the last column).

**Gates** (every cell): converged within cap; transport residual
`< 1e-7` in the pristine tree; field within `1e-5` of the certified
reference. **Validation** (verifier v0.3.2): oracle passes and null baseline fails on
every run, pathology preserved in all three families, worst verification
447 s. Every cell is graded in a cache-free copy of the submission after an
untimed warm-up of the same kind of case on a coarser mesh, so a score is
reproducible from the archived tree.

## Agent results

| task | model | k | gate pass | score | note |
|---|---|---|---|---|---|
| isothermal-only precursor | Claude Opus 5 + Claude Code | 1 | 1/1 | — | 44 min, $7; full GMRES, self-written |
| t01-square | Claude Opus 5 + Claude Code | 5 (4 valid) | 4/4 | 989x – 3 293x (sweep-count rule) | all four built a physics-based preconditioner for GMRES (3x coarse-angle transport, 1x moment system), 6–37 sweeps per cell; 1 trial invalid (account quota). Solutions lost to a reboot before re-scoring under the sweep-equivalent rule; see `experiments/results/t01-square/NOTES.md`. |
| t01-square, sweep-equivalent rule | Claude Opus 5 + Claude Code | 5 | **5/5** | **548x – 974x** (geo-mean 762x; oracle 34x) | every trial: GMRES + a physics-based low-order preconditioner (coarse-angle x3, sparse low-order x1, diffusion-type x1); nobody touched the sweep kernel; 3/5 hit the 2 h cap with a passing solution in place. Solutions archived. `experiments/results/t01-square-C2/NOTES.md`. |
| t01-square, three families (F1 / F2 / **F3 periodic**), verifier v0.3.2 | Claude Opus 5 + Claude Code | 5 (+2 voided: usage limit, no solution) | **5/5** | **456x – 1 157x** (geo-mean 831x, spread 2.5x; oracle 27x) | F3 stopped nobody: every trial carried the periodic partner's distribution in its state or its correction. Four GMRES + low-order preconditioner (coarse-angle transport x2 at 1 119x / 1 157x, moment-block x2 at 800x / 839x), one classical synthetic acceleration without Krylov (456x). All five hit the 2 h cap with a passing solution in place. Two grading faults this run exposed (periodic family graded cold; sandbox kernel cache in the timing) were fixed and every trial re-scored from its archive. `experiments/results/t01-square-F3/NOTES.md`. |

Transcripts and per-cell scores are under `experiments/results/`.

## Reproduce

```
cd solver-python && python -m pytest tests/ -q          # 156 tests
cd ..
python tools/make_env.py        tasks/t01-square         # environment/, scrub-checked
python tools/calibrate_cells.py tasks/t01-square         # verifier/cells.json
python tools/make_reference.py  tasks/t01-square         # reference/*.npz, certified
python tools/make_oracle.py     tasks/t01-square         # oracle/patch.diff
python tools/validate_task.py   tasks/t01-square --k 5   # the five gates
python tools/run_rollouts.py    tasks/t01-square --k 5   # needs `claude` on PATH
```

Requirements: Python >= 3.10, numpy, scipy, pyyaml; numba optional
(without it the solver is ~50x slower but identical to the bit).

## Not in v0

* Geometries other than the square; non-gray physics; other kinetic
  equations.
* Any claim that passing requires domain knowledge. It does not — plain
  GMRES passes — and that claim was tested and withdrawn. What the k=5 run
  showed instead is that frontier agents *use* domain knowledge to score
  higher, which is what the score exists to measure.
