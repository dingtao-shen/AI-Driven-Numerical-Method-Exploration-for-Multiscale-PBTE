# Make the solver converge across the whole stiffness range

`environment/` holds a working solver for the steady phonon Boltzmann
transport equation on a square: nodal discontinuous Galerkin in space,
discrete ordinates in angle, and an outer iteration that alternates a
transport sweep with a moment update.

Ten cases in `cases/` cover two boundary-condition families on the same
square and the same discretisation, each at five `(tau_R, tau_N)` pairs:

| family | walls | cases |
|---|---|---|
| `F1_*` | all four isothermal, `T = 1` on the north wall | five |
| `F2_*` | east and west diffusely reflecting (adiabatic), north hot, south cold | five |

The five pairs span the diffusive, transition, ballistic and hydrodynamic
regimes of the dual-relaxation-time model. The solver handles some of them
comfortably. On the stiff ones it slows down by orders of magnitude, and on
the stiffest it does not reach the answer at all inside its iteration cap.

**Your job: make every case converge — to the same answer the solver already
converges to where it can — and make it fast.**

## How it is graded

Grading has two parts: gates, which every case must pass, and a score.

**Gates** (all ten cases):

1. **It converges** inside its own `tmax`.
2. **The converged distribution solves the discrete transport system.** The
   grader computes the residual of that system directly from your converged
   distribution, with its own copy of the discretisation, and requires
   `||A f - b|| / ||b|| < 1e-7`. This does not go through your convergence
   test, so loosening a tolerance, damping the physics, or declaring
   convergence on some other quantity will not move it.
3. **The temperature field matches the stored reference** for that case to
   `1e-5` (relative, worst element), both as your run reports it and as
   recomputed from your converged distribution. The reference is the
   discrete fixed point of the shipped solver, certified independently. Note
   what this implies: stopping when successive iterates differ by `1e-8` is
   not the same as being within `1e-5` of the answer.

4. **The sweep is still there.** `solver.ctx.sweep(solver.mom, solver.vdf)`
   must remain callable and perform one full transport sweep from the
   current moments; the grader times it.

**Score**, only if every gate passes: the speed-up over the shipped solver
on the cases where it is slow, combined as a geometric mean per family and
then across families. Cost is measured in **sweep-equivalents**: the grader
times your `Case.from_yaml + Solver(case) + run()` end to end — setup
included — and divides by the time of **your own** transport sweep, timed by
the grader in the same single-threaded session (and capped by the grader's
own sweep time). Making the sweep itself faster therefore gains nothing —
the unit shrinks with the wall-clock — and making it slower is charged at
the grader's unit. What moves the score is work *outside* the sweep: a
coarse solve, a factorisation, a direct solve is charged at what it
actually costs, in sweep units. The shipped solver's cost is measured the
same way. Your reported iteration count and peak memory are shown alongside
the score but do not enter it.

The grader supplies its own copies of the case files and the mesh, so
changing the discretisation, the angular resolution or the stopping
tolerance in `cases/` will not affect grading.

## What the grader needs from your code

Keep these working — they are how the grader talks to your solver:

```python
from pybte import Case, Solver

case   = Case.from_yaml(path)      # must accept the shipped case files unchanged
solver = Solver(case)
record = solver.run()              # record.iterations, .converged, .temp
solver.vdf                         # the converged distribution, (ndir, n_tris, ndof)
```

`record.iterations` should count the transport sweeps your run performed; it
is reported next to the score. Anything else in the package is yours to
restructure — but note gate 4: the sweep stays callable as `solver.ctx.sweep`.

## Running things

```
cd environment
python -m pybte cases/F1_1e-2.yaml       # one case
python -m pytest tests/ -q               # the solver's own test suite
```

The test suite that ships with the solver pins the discretisation, the
quadrature, the sweep and the boundary treatment against known values, so
if a change of yours breaks the physics it will usually say so.

`README.md` in `environment/` describes the model, the two wall types and
the layout.

---

BENCHMARK-CANARY-a7f3c21e-DO-NOT-REPRODUCE-THIS-STRING
