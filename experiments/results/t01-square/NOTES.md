# t01-square — k=5, Claude Opus 5 + Claude Code, directory isolation

Run 2026-09-08. Grading at the time: gates + score with the score in
**fine-sweep count** (the rule in `task.md` the agents saw). The score rule
was changed to sweep-equivalents afterwards (see below); these trials could
not be re-scored under it.

| trial | gate | score (sweep count) | min | turns | cost | method (from the delivered code, read before it was lost) |
|---|---|---|---|---|---|---|
| 1 | PASS 10/10 | 3 293x | 108 | 105 | $18.97 | GMRES, **left-preconditioned by a coarse-angle transport operator**: the same DG discretisation on 16 ordinates, sparse LU 19 200 x 19 200 (1.76e6 nnz) factorised once. 6–8 sweeps on every cell. ~0.13 s per iteration against 1.1 ms per sweep; ~1 s per cell end to end; RSS 350 MB. |
| 2 | PASS 10/10 | 989x | 120 (cap) | — | — | GMRES **preconditioned by a discretised moment system** (energy and momentum, to reach the shear modes at small `tau_N`). 11–37 sweeps; ~0.35 s per cell; RSS 190 MB. Hit the 2 h cap; transcript lost; graded as left. |
| 3 | PASS 10/10 | 2 605x | 61 | 77 | $11.08 | GMRES **preconditioned by a coarse-angle transport operator** (independent implementation of trial 1's idea). 8–11 sweeps; ~0.43 s per cell; RSS 300 MB. |
| 4 | PASS 10/10 | 2 271x | 91 | 100 | $17.00 | GMRES **preconditioned by a coarse-angle transport operator** (8 equally spaced in-plane directions, weights and speed chosen to match the first two moments). 9–13 sweeps; ~0.5 s per cell; RSS 330 MB. |
| 5 | — | — | 17 | 23 | $3.15 | **Invalid.** The account's session limit was reached ("You've hit your session limit"); the sandbox was graded as the unmodified environment. Environment friction, excluded from the count. |

**Valid trials: 4. Gate pass rate: 4/4. Scores 989x – 3 293x (spread 3.3x).**

Reference points, same rule: source iteration 33 .. > 200 000 sweeps; the
Krylov oracle (unpreconditioned full GMRES) 21 .. 1 104 sweeps, score 42x.

## What the trials showed

* Every valid trial went beyond plain GMRES to a **physics-based
  preconditioner** — three chose a coarse-angle transport operator, one a
  moment system — and every one landed on an iteration count nearly
  independent of the Knudsen number. That is the property the v0.1 proposal
  said no agent would reach without domain knowledge. They reached it with
  domain knowledge, in 1–2 hours each.
* The preconditioners differ in cost, and the sweep-count score is blind to
  it: trial 1 spends ~120 sweeps' worth of time in each low-order solve,
  trial 2 ~30. By sweep count trial 1 ranks first; by wall-clock per cell
  trial 2 does. This is why the score rule was changed.
* The 2 h agent budget is binding (trial 2). The account quota is a
  failure mode of the harness, not of the task (trial 5).

## Score rule change and what it did to the record

The score is now in **sweep-equivalents**: end-to-end wall-clock of
`Case + Solver + run` divided by the time of one fine sweep, both timed by
the verifier in-session; the unaccelerated solver's cost is timed the same
way. Under that rule the oracle scores **19.4x** (F1 >= 23.5x, F2 >= 16.0x)
instead of 42x, because its setup now counts.

The five sandboxes were held in the session scratch directory and were
lost when the machine rebooted on 2026-09-09 (14:25). They cannot be
re-scored, and the code the agents wrote survives only in the descriptions
above. `tools/run_rollouts.py` now archives every trial's solution — a
unified diff and a tarball of the environment — into this directory at
grading time. A fresh run under the sweep-equivalent rule is required for a
publishable number; the four scores above are the sweep-count rule's.
