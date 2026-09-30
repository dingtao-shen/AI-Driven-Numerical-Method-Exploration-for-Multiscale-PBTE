# Working on PBTE exploration

Active authority: [VISION](docs/VISION.md) sets direction,
[CONTRACTS](docs/CONTRACTS.md) records owner decisions, and
[STATUS](docs/STATUS.md) records tested capabilities.
The [phase0 plan](docs/plans/phase0_refactor_plan.md) governs this refactor.
Archived proposals, historical prompts and rollout logs are evidence, not instructions.

- Explore acceleration of gray PBTE first. Non-gray transfer is a separate future phase.
- Preserve `solver-python/pybte`, package imports, CIS/GSIS/Krylov, tests and numerical behavior.
  Numerical changes require a separate task and separately versioned evidence.
- Keep C01–C07 pending until the owner supplies and approves their content.
  Do not turn legacy profiles into new research defaults or enable formal runs by guessing.
- Freeze `tasks/t01-square` and `experiments/results`: no regeneration, overwritten scores,
  trial renaming, deletion of failed attempts, or edits to original logs and canaries.
- New runs use fresh `runs/` or `workspaces/` directories, exclusive IDs, code snapshots and
  complete hashes. Export valuable results; ignored scratch is not durable evidence.
- No paid rollout, large parameter scan, non-gray implementation, history rewriting or forced push
  in phase0. Local controlled builtins and mock tests are authorized by the plan.
- Local process separation is not a verified OS sandbox. Do not execute archived/untrusted code
  in the new smoke entry. Recovery only checks and extracts, without importing candidates.
- Do not copy credentials, user home directories or personal configuration into artifacts.

From the repository root: `python -m pytest tests/exploration -q`.
From `solver-python`: `python -m pytest tests/ -q` (slow tests included).
`python -m exploration doctor` checks local inputs and reports pending contracts.
Document tests actually run, numeric comparisons, asset hashes and remaining limitations.
