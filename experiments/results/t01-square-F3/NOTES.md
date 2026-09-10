# t01-square v0.3 (three families) — k=5, Claude Opus 5 + Claude Code

Run 2026-09-10 on the three-family task (`a711541`): F1 isothermal, F2
diffusely reflecting sides, F3 periodic sides; five `(Kn_R, Kn_N)` points
each; cap 6 000; own-sweep-unit score; grading pinned to one thread. Oracle
under this rule: **28.3x** (F1 >= 35.1x, F2 >= 34.9x, F3 >= 18.4x; `oracle_score.json`).
F3 is its weakest family too: the periodic-trace block makes its Krylov vectors
seven times longer (27 600 entries) and the basis 617 MB against ~220 MB elsewhere. Sandboxes under
`~/rollout-sandboxes/t01-square-F3/`; solutions archived here as diff +
tarball; transcripts as stream-json (`trial_XX_transcript.jsonl`), which
survive a timeout.

The question this run asks: **does the periodic family stop an agent whose
first instinct is Krylov on the moments?** In the first k=5, the unaided
GMRES solution declared periodic faces inapplicable and fell back to source
iteration, failing the stiff cells.

| trial | gate | score | F1 / F2 / F3 | min | method (from the archived diff) |
|---|---|---|---|---|---|
| 1 | PASS 15/15 | 689x | 747x / 859x / 510x | 120 (cap; solution in place) | GMRES on a state of moments **+ the distribution behind each periodic face + wall emission**, preconditioned by low-order (DSA/MSA-type) operators on the moment block. F3 costs ~2 sweep-equivalents per sweep against ~1.2 on F1: the periodic-trace part of the state is carried but not preconditioned. |
| 2 | | | | | |
| 3 | | | | | |
| 4 | | | | | |
| 5 | | | | | |

## Observations so far

* **The periodic family did not stop trial 1.** The agent identified exactly
  what the oracle does — that a periodic face's partner is part of the state
  the sweep reads — and carried it in the Krylov vector. F3 is the costliest
  family for it (510x against 747x / 859x), which the per-family score shows,
  but it passes.
* Stream-json kept the transcript through the timeout (2 207 events).
