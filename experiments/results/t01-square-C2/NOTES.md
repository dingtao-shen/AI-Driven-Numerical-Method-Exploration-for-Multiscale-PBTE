# t01-square — k=5 under the own-sweep-unit rule, Claude Opus 5 + Claude Code

Run 2026-09-10. Score rule (`task.md`, commit `a0765ac`): speed-up in
sweep-equivalents where the unit is **min(the submission's own fine sweep,
the verifier's fine sweep)**, both timed single-threaded in-session; the
submission's cost is `Case + Solver + run` end to end, setup included. A
faster sweep gains nothing, a slower one is charged in the verifier's unit,
only work outside the sweep moves the score. Gate 4 requires
`solver.ctx.sweep` to stay callable. Reference: the Krylov oracle
(unpreconditioned full GMRES) scores **34.0x** under this rule.

| trial | gate | score | min | turns | cost | method (from the archived diff) | sweeps / sweep-eq per cell |
|---|---|---|---|---|---|---|---|
| 1 | PASS 10/10 | 548x | 120 (cap) | — | — | GMRES + **coarse-angle low-order transport operator** (own velocity mesh, own sweep order). Solution was in place when the cap hit; transcript lost. | 11–21 / 22–108 (the two hydrodynamic cells cost 85–108: its low-order solve is expensive there) |
| 2 | PASS 10/10 | 921x | 120 (cap) | — | — | GMRES + **sparse low-order operator** (scipy LU) on the macroscopic unknowns. Leanest overhead of the four. Solution in place at the cap; transcript lost. | 14–25 / 20–37 |
| 3 | PASS 10/10 | 974x | 84 | 87 | $15.64 | **Right-preconditioned GMRES** + coarse-angle low-order operator (`krylov.py`, `loworder.py`). | 11–26 / 18–38 |
| 4 | PASS 10/10 | 862x | 110 | 111 | $21.96 | GMRES + **coarse-angle low-order operator** (`coarse_directions`, `LowOrderOperator`). | 13–24 / 21–39 |
| 5 | | | | | | | |

## What changed with the rule

* **Nobody touched the sweep kernel.** Own-sweep times 8.4–11 ms against the
  verifier's 8.6–10.6 ms on every cell of every trial. The C-kernel detour of
  the previous run did not recur; there was no longer anything to gain from it.
* **All four built the same kind of thing as under the sweep-count rule** — a
  physics-based preconditioner for GMRES — and the scores now carry its cost:
  548x – 974x here against 989x – 3 293x by raw sweep count, because the
  low-order solve is charged at 1.3–2.4 sweep-equivalents per sweep (up to 8
  on trial 1's hydrodynamic cells).
* **Spread 1.8x** across four trials of one model (974 / 548). The score
  separates the same idea implemented with more or less overhead, which is
  what it was changed to do; it ranks trial 2's lean sparse solve above
  trial 1's heavier one although trial 1 uses fewer sweeps.
* **Against the oracle**: 16x – 29x better than unpreconditioned GMRES under
  the same rule. That gap is the value of the preconditioner, measured at its
  real cost.
* **The 2 h cap** was hit by two of four and neither was hurt: both had
  working solutions in place and were polishing. Transcripts are lost on
  timeout (`--output-format json`); `stream-json` would keep them.
