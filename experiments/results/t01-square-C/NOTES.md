# t01-square — k=5 under the sweep-equivalent rule, Claude Opus 5 + Claude Code

Run started 2026-09-09. Score = speed-up in sweep-equivalents (end-to-end
wall-clock of `Case + Solver + run` over the time of one fine sweep, both
timed by the verifier in-session, setup included), as stated in `task.md`.
Sandboxes under `~/rollout-sandboxes/t01-square-C/`; each trial's solution
is archived here as a diff and an environment tarball at grading time.

Grading was pinned to one thread (`OMP/MKL/OPENBLAS/NUMBA_NUM_THREADS=1`)
from trial 2 on; trial 1 produced no runnable tree, so no score in this run
was taken under unpinned grading.

| trial | gate | score | min | turns | cost | what the delivered tree contains |
|---|---|---|---|---|---|---|
| 1 | FAIL 0/10 | — | 120 (cap) | — | — | **Timed out mid-edit; every cell crashed** (`_bte.so: undefined symbol: bte_set_threads`). The diff shows compiled **C kernels** for the sweep (`_csrc/bte.c`, a ctypes binding with thread control, a build script), plus `krylov.py` and a low-resolution preconditioner `lowres.py`. The `.so` was stale against the C source when the cap hit, and the binding's "fall back if the library will not load" guard did not cover a missing symbol. Graded as left. |
| 2 | | | | | | |
| 3 | | | | | | |
| 4 | | | | | | |
| 5 | | | | | | |

Reference: the Krylov oracle under this rule scores **33.8x** with grading pinned to
one thread (F1 >= 33.2x, F2 >= 34.4x; `oracle_score.json`), 19.4x unpinned -- the
single-threaded sweep is a larger unit, so the oracle's non-sweep overhead weighs less.

## Observations so far

* **The rule changed what the agent optimised.** Under the sweep-count rule
  every valid trial built a preconditioner; under the wall-clock-based rule
  the first trial spent its budget compiling the sweep and threading it —
  legitimate under the letter of the task, off the point of the benchmark,
  and fatal when the cap arrived mid-build. A score that pays for
  non-sweep work at its true cost also pays for a faster sweep, and a faster
  sweep is easier to get from a compiler than from a better iteration.
* The verifier now pins grading to one thread so that a parallel sweep does
  not score core count. It does not, and cannot, stop a submission from
  making the single-threaded sweep faster; whether that should count is a
  design question left open here.
* Transcripts are lost on timeout (`--output-format json` emits only at the
  end); `stream-json` would keep them.
