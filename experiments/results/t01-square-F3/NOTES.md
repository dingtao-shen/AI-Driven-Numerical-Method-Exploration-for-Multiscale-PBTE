# t01-square v0.3 (three families) — k=5, Claude Opus 5 + Claude Code

Run 2026-09-10 → 09-14 on the three-family task: F1 isothermal, F2
diffusely reflecting sides, F3 periodic sides; five `(Kn_R, Kn_N)` points
each; cap 6 000; own-sweep-unit score; grading pinned to one thread; 2 h
agent budget. **Every trial and the oracle were graded uniformly under
verifier v0.3.2 (`ac1d315`) from the archived tarballs, on an idle machine,
2026-09-14** (`tools/rescore_rollouts.py`). Oracle under this rule:
**27.0x** (F1 34.8x, F2 32.4x, F3 17.5x; `oracle_score.json`). Sandboxes
under `~/rollout-sandboxes/t01-square-F3/`; solutions archived here as diff +
tarball; transcripts as stream-json.

The question this run asks: **does the periodic family stop an agent whose
first instinct is Krylov on the moments?** In the first k=5, the unaided
GMRES solution declared periodic faces inapplicable and fell back to source
iteration, failing the stiff cells.

## Result: 5/5 pass, 456x – 1 157x (geometric mean 831x, spread 2.5x)

| trial | gate | score | F1 / F2 / F3 | sweeps per cell | non-sweep cost (sweep-eq per cell) | method (from the archived diff) |
|---|---|---|---|---|---|---|
| 1 | PASS 15/15 | 800x | 782 / 890 / 735 | 12 – 54 | 4 – 15 | GMRES on a state of moments **+ the distribution behind each periodic face + wall emission**, preconditioned by low-order (DSA/MSA-type) operators on the moment block. |
| 2 | PASS 15/15 | 839x | 835 / 905 / 782 | 15 – 30 | 3 – 20 | Same structure ("the sweep ordering treats a periodic face as a boundary, so that inflow is lagged"); GMRES preconditioned by a coarse low-order operator (`loworder.py`). |
| 3 | PASS 15/15 | 456x | 434 / 516 / 423 | 8 – 27 | 27 – 47 | **Not Krylov.** Classical synthetic acceleration: the same DG transport operator with nine directions, factorised once, applied as a correction after every sweep; the correction is also injected into the distribution (P1 form) at the elements a sweep reads stale — wall emission and periodic inflow — "with them lagging the iteration diverges outright on F2 and F3". Stops on the per-element size of the correction. Fewest sweeps, highest overhead. |
| 4 | PASS 15/15 | 1 157x | 1 138 / 1 296 / 1 050 | 7 – 20 | 7 – 24 | GMRES on the full state; preconditioner = Galerkin projection of the DG transport operator onto angular polynomials of degree 1 (degree 2 when `tau_N << tau_R`, "because the damping of the momentum modes is set by the fourth angular moments"), sparse LU once; converged distribution assembled by linearity, no closing sweep. |
| 5 | PASS 15/15 | 1 119x | 1 030 / 1 209 / 1 126 | 7 – 20 | 8 – 17 | Right-preconditioned GMRES on the full state; preconditioner = coarse-angle synthetic transport solve (same DG discretisation and boundary treatments, a handful of directions, scattering / wall / periodic coupling folded in), SuperLU once. Added a wall-temperature override to `sweep()` to evaluate the homogeneous part of the affine map. |

All five hit the 2 h cap with a passing solution in place (stream-json kept
every transcript). Peak RSS 252 – 317 MB.

**Voided, not counted:** the original trials 4 and 5 (2026-09-11) were ended
by the account's usage limit (API 429) after 16 min / 28 turns and 4 s /
1 turn, with no solution written; their meta, transcript and (null) score are
under `voided/`. `run_rollouts.py` now marks such trials voided and
`--trials 4 5` replaced them on 2026-09-14.

## Observations

* **The periodic family did not stop anyone.** All five identified what the
  oracle does — that a periodic face's partner is part of the state the
  sweep reads — and either carried it in the Krylov vector (1, 2, 4, 5) or
  wrote the correction into it (3). F3 is not systematically the weakest
  family for the agents (it is for the oracle, whose periodic-trace block is
  unpreconditioned).
* **The score separates methods that all pass.** Coarse-angle-transport
  preconditioners for GMRES (4, 5): 1 119x – 1 157x. Moment-block
  preconditioners for GMRES (1, 2): 800x – 839x. Synthetic acceleration
  without Krylov (3): 456x — fewest sweeps (8 – 27) but 27 – 47
  sweep-equivalents of low-order work per cell. The ordering is by sweep
  count and non-sweep cost together, which is what the unit was designed to
  charge. Nobody optimised the sweep kernel.
* **Noise floor of one grading: about ±3%.** Trials 4 and 5 were graded at
  the end of their runs and again 30 min later from the archive: 1 148.8x →
  1 157.0x and 1 157.4x → 1 119.4x. The machine's sweep time also drifted
  ~10% between days (8.9 → 9.8 ms); the ratio cancels most of it. The 2.5x
  spread across trials is well above this.

## Two grading faults this run exposed (fixed before the numbers above)

Both were found by re-scoring archived solutions, which is what the archive
is for. The original in-sandbox scores (verifier v0.3.0) were trial 1
689x (747 / 859 / 510), trial 2 740x (810 / 927 / 540), trial 3 424x
(424 / 530 / 339).

1. **The periodic family was graded cold** (v0.3.0). The warm-up solve
   switched every wall to thermalising but kept the `Master`/`Slave` names
   while dropping their offsets; the mesh reader pairs periodic faces by
   those names and threw `unpaired periodic faces` on every F3 case, for the
   oracle and every agent. The graded F3 run then paid the JIT: ~0.1 – 0.3 s
   on runs of 0.2 – 0.3 s. The "F3 costs ~2 sweep-equivalents per sweep"
   reading of the first two trials was this, not the periodic state.
2. **The sandbox's compiled-kernel cache was part of the score** (v0.3.1).
   pybte's kernels are disk-cached; graded in the sandbox it was developed
   in, a submission had every kernel pre-compiled, while the same tree
   re-scored from its archive compiled the reflecting-wall kernels inside
   the first F2 cell's timing (0.18 s → 1.2 s; trial 1 F2 859x → 483x).

v0.3.2: the warm-up runs the graded case's own wall types on the coarse mesh
of the same family (offsets kept), the verifier refuses to score a cell whose
warm-up did not run, and every cell runs in a copy of the submission stripped
of `__pycache__`. Validation under v0.3.2: oracle PASS 2/2 (142 s), null
FAIL 2/2 (447 s), pathology preserved in all three families. Trial 1's
archive now scores the same whether graded in place or from the tarball.
