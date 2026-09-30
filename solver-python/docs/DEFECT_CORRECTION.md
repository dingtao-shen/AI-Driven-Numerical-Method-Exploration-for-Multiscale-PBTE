# Experimental correction of the observed GSIS fixed-point discrepancy

**Status: research code, off by default.** Turning it on changes every GSIS
iteration count. Nothing in the shipped cases, the regression suite, or the
validation numbers uses it.

Read `LIMITATIONS.md` #1 first. Tables are historical measurements of this
implementation; causal explanations below are hypotheses requiring further
verification, not universal conclusions about GSIS.

## What is being repaired

The last step of a GSIS iteration blends the synthetic macroscopic solution
with the moment of the kinetic solution:

```
M^{n+1} = beta M* + (1 - beta) M(f)          f = f^{n+1/2}
```

Write the bracket it blends in as

```
o = M* - M(f)
```

`M*` comes from the HDG synthetic system, `M(f)` from the upwind-DG kinetic
solution. If the two discretisations were exact moments of one another `o`
would vanish at convergence and the blend would be inert there. The recorded implementation measured a stationary nonzero offset
`|o_T| = 1.49e-01` at `tau_R = 1e-1`. Establishing its cause requires
checking the discrete equations and implementation, not just their names.

## The repair

Carry an estimate `d` of the stationary offset and subtract it before the
blend:

```
M^{n+1} = M(f) + beta (o - d)
d       <- d + omega (o - d)
```

At a fixed point `d = o`, the blend vanishes, and `M = M(f)` — **exactly
CIS's fixed-point condition, for any `beta` and any macroscopic
discretisation.** `omega = 0` leaves `d` identically zero and reproduces
the shipped scheme bit for bit.

Everything from the sweep through the HDG solve to the blend is linear in the
state, so `d -> o(d)` is an **affine** map. That is what licenses Anderson
acceleration on it below, under the assumed linear state representation; this does not establish
the full causal diagnosis or robustness beyond the tested cases.

## Update policy is the whole story

### Refreshing every iteration destroys the acceleration

`scheme.defect_every: 1`, `tau_R = 1e-1` (CIS: 557 iterations, mass
`0.2499301936`):

| `defect_omega` | iterations | gap vs CIS | transport residual | mass |
|---|---|---|---|---|
| 0 (published) | **49** | 1.687e-02 | 3.16e-03 | 0.2499964771 |
| 0.1 | 4297 | 1.7e-11 | 4.18e-13 | 0.2499301936 |
| 0.3 | 1869 | 5.6e-12 | 1.65e-13 | 0.2499301936 |
| 0.5 | 1354 | 3.3e-12 | 1.15e-13 | 0.2499301936 |
| 1.0 | diverges (overflow) | — | — | — |

The fixed point is repaired perfectly — the mass agrees with CIS to all ten
printed digits and the transport residual falls to round-off — and the
acceleration is gone: 1354 iterations against CIS's 557.

The reason is the lag. Linearise about the fixed point with `G` the CIS map
and `H` the synthetic map:

```
e^{n+1} = [Gamma + beta(Eta - Gamma)] e^n  -  beta(Eta - Gamma) e^{n-1}
```

The acceleration term has become a *difference of successive offsets*, which
is second-order small near convergence, so the iteration degenerates to CIS
plus a destabilising lag. With `omega = 1` the companion matrix leaves the
unit disc and it diverges outright.

### Freezing the defect leaves the iteration matrix alone

`scheme.defect_every: 0` refreshes `d` only once the inner iteration has
converged with `d` held fixed. Between refreshes the iteration matrix is
**exactly** the shipped scheme's, so the inner convergence rate is
unchanged; only the fixed point moves. The cost moves into an outer loop.

Gap after a bounded number of outer cycles, `tau_R = 1e-1`, `omega = 1`:

| outer cycles | iterations | gap vs CIS | transport residual | x CIS cost |
|---|---|---|---|---|
| 0 | 49 | 1.687e-02 | 3.16e-03 | 0.09 |
| 1 | 85 | 1.215e-02 | 9.80e-04 | 0.15 |
| 3 | 154 | 8.618e-03 | 4.36e-04 | 0.28 |
| 8 | 317 | 3.862e-03 | 1.76e-04 | 0.57 |
| 20 | 680 | 6.665e-04 | 3.29e-05 | 1.22 |
| 40 | 1223 | 4.965e-05 | 2.20e-06 | 2.20 |

The outer loop contracts at `rho = 0.950` (`omega = 1`) or `0.973`
(`omega = 0.5`) — against CIS's own `rho ~ 0.977` at this Knudsen number.
The similar contraction factors suggest a slow mode shared with CIS in this
case. Identifying that mode and attributing it to a specific discretization
requires further analysis; the measured cost and convergence table are retained.

### Anderson on the outer loop

The outer map is affine, so Anderson (type II, `scheme.defect_anderson` =
window) is exact arithmetic on it rather than a heuristic. At `tau_R = 1e-1`,
`omega = 1`, to full convergence:

| outer policy | outer cycles | iterations | gap vs CIS | transport residual |
|---|---|---|---|---|
| plain relaxation | 317 | 4480 | 3.96e-12 | 4.42e-14 |
| Anderson, window 8 | **69** | **899** | 2.69e-12 | 1.78e-14 |

Five times fewer outer cycles for a few least-squares solves on a
`3*NDOF_TRI x N_TRIS` vector, which is free next to the sweeps.

## Across Knudsen number

`tau_R = 1e-2`, 200 elements, `DEG = 3`, `tol = 1e-13`. The reference is CIS
run to round-off; `gap` is `||T - T_CIS||_inf / ||T_CIS||_inf`.

| scheme | iterations | outer | transport residual | gap vs CIS | `int T dA` | wall |
|---|---|---|---|---|---|---|
| CIS | 34 846 | — | 9.55e-14 | 0 | 0.2499935255 | 337 s |
| GSIS (published) | **43** | — | 2.28e-04 | 3.72e-03 | 0.2499619359 | 1 s |
| GSIS + defect | 60 000 (not conv.) | 3 708 | 5.31e-09 | 6.73e-06 | 0.2499925063 | 761 s |
| GSIS + defect + AA(8) | 5 317 | 417 | **1.13e-13** | **1.05e-10** | 0.2499935255 | 91 s |
| GSIS + defect + AA(16) | **4 616** | 383 | 2.26e-13 | 5.03e-10 | 0.2499935256 | **78 s** |

**This is the case for the method.** At `tau_R = 1e-2` the repaired scheme
reaches the exact discrete kinetic fixed point — mass agreeing with CIS to
ten digits, transport residual at round-off — in **7.5x fewer iterations and
4.3x less wall time than CIS**.

The advantage is a function of Knudsen number, because the two costs scale
differently:

| `tau_R` | CIS | GSIS + defect + AA(8) | ratio |
|---|---|---|---|
| 1e-1 | 557 | 899 | 0.6x (worse) |
| 1e-2 | 34 846 | 5 317 | **6.6x** |

CIS grows roughly as `tau_R^-1.8` here, the repaired scheme as `tau_R^-0.7`.
It is **not** Knudsen-independent the way the shipped scheme is (49 and 43
iterations at `1e-1` and `1e-2`) — the outer loop still has to converge the
slow diffusive mode that the offset carries. It merely converges it much
faster than source iteration does.

### Deep in the diffusive regime the comparison inverts

`tau_R = 1e-3`, everything capped at 60 000 iterations. **CIS never gets a
usable answer**, so there is nothing to take a gap against; the transport residual remains a useful observation. A small residual alone
does not certify field accuracy for an ill-conditioned system; reference
certification and field checks must be stated separately.

| scheme | iterations | outer | transport residual | `int T dA` | wall |
|---|---|---|---|---|---|
| CIS | 60 000 (not conv.) | — | 7.63e-06 | 0.1336293140 (89% wrong) | 578 s |
| GSIS (published) | **93** | — | 5.81e-06 | 0.2499975333 | 1 s |
| GSIS + defect | 60 000 (not conv.) | 2 213 | 1.38e-08 | 0.2499981845 | 758 s |
| GSIS + defect + AA(8) | 60 000 (not conv.) | 3 573 | **4.90e-11** | 0.2499991791 | 1010 s |

Two things to read off. First, the published GSIS at 93 iterations already
has a *smaller* transport residual than CIS has after 60 000 — the
displacement shrinks with `tau_R`, so there is progressively less to repair
just as CIS becomes progressively less able to serve as a reference. Second,
the repaired scheme still buys five orders of magnitude on the residual over
the uncorrected implementation. That residual improvement is not itself a
complete field-accuracy certification.

## Evidence and limits

The retained tables show that the experimental correction can reduce the observed
fixed-point discrepancy in sampled cases, sometimes at substantially higher
iteration cost. The `tau_R = 1e-2` comparison reported 4.3x lower wall time than
its CIS comparison; those historical timings are not a new-protocol score.
The correction is off by default. Its diagnosis, broader parameter robustness,
and compatibility with alternative synthetic discretizations remain unverified.
It is not automatically the preferred reference generator; certification must
be recorded case by case, and the existing Krylov path is documented separately.

## Configuration

```yaml
scheme:
  accflag: 1
  defect_omega: 1.0     # 0 = the shipped blend; relaxation onto the offset
  defect_every: 0       # 0 = refresh only once the inner iteration converges
                        # N = refresh every N iterations (N=1 kills the speed-up)
  defect_anderson: 8    # Anderson window on the outer sequence; 0 = plain
```

Historical experimental configuration: `defect_omega: 1.0`, `defect_every: 0`,
`defect_anderson: 8`. Do not use `defect_every: 1` — see the first table.

Diagnostics land in `RunRecord.diagnostics`: `defect_omega`, `defect_every`,
`defect_anderson`, `defect_outer` (outer cycles taken) and `defect_norm`
(`max|d|`). `Solver.acc.defect_steps` is the outer convergence history.

Reproduce with `tools/defect_study.py`:

```
python tools/defect_study.py --tau-r 1e-2 --omega 1.0 --every 0 --anderson 8
```
