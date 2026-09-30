# Known limitations

What this solver cannot do, or does in a way you need to know about. Each
entry is a property of `pybte` **now** — none is a "to do".

Defects that were found in the Fortran original and *fixed* during the port
are not listed here; they no longer affect anything. They are recorded in the
git history at the `Port ACC_2D2V_LinearCallawayModel` commit and the ones
following it, together with the reference sources they were found in.

---

## 1. Recorded fixed-point differences in the current GSIS/CIS implementation

**Evidence scope:** the tables describe the shipped implementation and tested cases.
Causal interpretation and discrete compatibility require further verification;
these observations do not prove a universal GSIS property or exclude a bug.

The last step of a GSIS iteration blends the synthetic macroscopic solution
`M*` with the moment `M(f)` of the kinetic solution, under a local damping
factor:

```
M^{n+1} = beta M*  +  (1 - beta) M(f)
beta    = min(tau_R/h_l, tau_thr) / (tau_R/h_l)
```

At a fixed point the blend requires `beta (M* - M(f)) = 0` to satisfy
CIS's moment consistency condition. The observed discrepancy motivates checking
the assembled macro/kinetic systems and implementation. Different discretization
names alone do not prove that this bracket must be nonzero. The damping factor
and the recorded ballistic observations are retained (see #5).

| | CIS | GSIS |
|---|---|---|
| iterations (`tau_R = 1e-1`) | 557 | 49 |
| transport residual at the answer | `4.3e-14` | `3.2e-03` |
| `int T dA` | 0.249930 | 0.249996 |

`||T_CIS - T_GSIS||_inf / ||T_CIS||_inf = 1.7e-2`. It is non-zero at **every**
Knudsen number tested, peaking around `tau_R = 1e-1`.

### Recorded dependence on beta

`tau_R = 1e-1`, 200 elements, `DEG = 3`:

| `flow.tau_thr` | `beta` range | GSIS iterations | gap vs CIS | transport residual |
|---|---|---|---|---|
| 0.01 | 0.0006–0.013 | 470 | **1.0e-03** | 6.1e-05 |
| 0.1 | 0.006–0.13 | 203 | 6.5e-03 | 4.9e-04 |
| **1.0** (default) | 0.06–1.0 | **49** | 1.69e-02 | 3.2e-03 |
| 10 | 0.61–1.0 | 57 | 1.97e-02 | 7.9e-03 |
| 1000 (`beta == 1`) | 1.0 | 91 | 1.98e-02 | 8.5e-03 |

In this sampled damping sweep, smaller beta reduced the gap and increased
iteration cost. This does not exclude other compatible accelerated formulations.

### Recorded refinement study

`tools/fixed_point_study.py`:

| varying | values | gap |
|---|---|---|
| `DEG` = 1, 2, 3 | fixed `tau_R = 1e-1`, 200 elements | 1.60e-2, 1.65e-2, 1.69e-2 |
| elements = 50, 200, 800 | fixed `tau_R = 1e-1`, `DEG = 3` | 1.66e-2, 1.69e-2, 1.43e-2 |
| `tau_R` = 1, 1e-1, 1e-2 | fixed discretisation | 1.15e-2, 1.69e-2, 3.72e-3 |

The sampled 16-fold element increase reduced the gap by about 14%.
These finite observations do not establish an asymptotic refinement result or
show that changing tau_R is the only way to achieve consistency.

### What follows

The recorded implementation cannot be treated as a reference for CIS's discrete
answer without checking the target problem. Correctness should use independently
computed observables under an explicitly stated protocol:

* the transport residual (`iteration.true_residual: true`) — necessary but,
  at small `Kn`, not sufficient on its own: the operator is ill-conditioned
  there, and a residual of `4e-12` was measured alongside a `1.5e-5` error in
  `int T dA`. Pair it with a field comparison against a certified reference;
* the analytic Fourier limit at small `Kn` (`pybte.analytic`);
* the iteration counts, which are unaffected and are the real phenomenon —
  16 836 against 27 at `tau_R = 1e-2`.

### The accelerated path that keeps the fixed point

`scheme.method: krylov` solves the source-iteration system `(I - T) u = g`
by GMRES, one sweep per product, and converges to **source iteration's own
fixed point** to round-off: transport residual `1e-15`, agreement with CIS
to `4e-11 .. 1.6e-7` wherever CIS can reach the answer. On the shipped
square at `tau_R = 1e-2` it takes 228 sweeps against CIS's 16 830; at
`1e-3`, 992 where CIS does not converge in 200 000. Its iteration count grows across the sampled Kn values. Periodic faces are
supported through the extended periodic partner state; unrepresented broken
sweep cycles still raise. See `KRYLOV.md` and the Krylov tests.

`scheme.defect_omega` is an earlier, experimental repair of GSIS itself
(`docs/DEFECT_CORRECTION.md`); it also lands on CIS's fixed point but at
3–8x the Krylov path's cost, and stays off by default.

## 2. Do not trust a CIS residual

The stopping criterion measures the **step** between iterates, not the
**error**. Source iteration converges linearly, `r_{n+1} = rho r_n`, so the
error remaining when you stop is the sum of all future steps,
`r rho/(1 - rho)` — and `rho -> 1` as the medium becomes optically thick.

| `tau_R` | CIS iterations | reported residual | absolute error in `int T dA` (`~0.25`) |
|---|---|---|---|
| 1e-2 | 16 836 (converged) | 1.0e-08 | 1.0e-05 |
| 1e-3 | 200 000 (truncated) | 1.5e-06 | **4.4e-02** |
| 1e-4 | 200 000 (truncated) | 2.5e-06 | **2.2e-01**, i.e. 89% wrong |

Directly measured at `tau_R = 1e-2` (200 elements, `DEG = 3`, `20 x 40`),
stopping at `tol = 1e-4`: 2 599 iterations, reported residual `9.99e-05`.
The reference is the **same scheme on the same discretisation** converged to
`tol = 1e-12` (31 245 iterations) — not the analytic limit, which would mix
discretisation error into an iteration measurement.

| quantity | relative error |
|---|---|
| `int T dA` | `1.28e-01` |
| cell-average `T`, area-weighted L2 | `1.08e-01` |
| cell-average `T`, L-infinity | `7.37e-02` |

Each is **three orders of magnitude** above the reported residual. The error
sits in the bulk — worst element at `(0.56, 0.62)`, corner elements twenty
times better — which is what an unconverged diffusive mode looks like. It is
not a boundary artefact, and the converged solution has no overshoot: `T`
stays inside `[0, 0.9884]` on every DOF.

`RunRecord.error_estimate` fits `rho` to the tail of the residual history and
reports `r rho/(1-rho)`; it lands within an order of magnitude of the truth.
`RunRecord.contraction` gives `rho` itself. **Look at these before believing a
CIS result.**

Note that `iteration.true_residual: true` does *not* expose this. For CIS it
measures the same moment change the iterate residual does (`1.52e-06` against
`1.52e-06` on the truncated `tau_R = 1e-3` run). What it is for is #1.

The historical stopping rule is retained for compatibility; new correctness
and same-accuracy baseline policies await C03/C04.

## 3. Acceleration variant B diverges

`scheme.acc_variant: B` rescales the momentum and stress equations by `TAU_R`.
It is algebraically the same system as variant A and was expected to be better
conditioned as `TAU_R -> 0`. It diverges geometrically from the first
iteration at every Knudsen number tested — `mass` goes 5.8, -193, 6772,
-2.4e5, reaching NaN within about 30 iterations — and the divergence is not a
conditioning artefact, since the two global matrices have comparable
1-norm condition estimates:

| `tau_R` | cond₁ A | cond₁ B | A iterations | B |
|---|---|---|---|---|
| 1e+00 | 2.5e7 | 4.9e7 | 22 | diverges |
| 1e-01 | 2.8e5 | 6.2e4 | 30 | diverges |
| 1e-02 | 1.6e4 | 1.4e5 | 27 | diverges |
| 1e-03 | 3.2e4 | 4.0e6 | 50 | diverges |

This is the behaviour of the variant as it was written in the original solver,
reproduced rather than worked around. It was never the variant that solver's
build system compiled. **Use variant A**, which is the default; B is kept
because the comparison is scientifically interesting and because a future
investigation may want the starting point.

## 4. Angular quadrature is accurate to ~1e-13, not to machine precision

The Gauss–Legendre routine is a transcription of the original's, which forms
the Legendre derivative as

```
Lp = (n+1) * (P_{n-1} - y P_n) / (1 - y^2)
```

where the identity calls for `n`, not `n+1`. The weight formula carries the
same `(n+1)^2/n^2` factor, so the weights stay consistent with whatever roots
come out — but the Newton step is scaled by `n/(n+1)`, which turns quadratic
convergence into linear convergence with ratio `1/(n+1)`. With a step-size
exit test (`|dy| < 1e-13`) the roots stop about `1e-13/n` short.

Measured: nodes differ from `numpy.polynomial.legendre.leggauss` by up to
`1.8e-14` (n=2) to `4.4e-16` (n=20); `sum(w)` on `[-1,1]` is short of 2 by
`~2e-13`; `sum(DOMEGA)` misses `4 pi` by `3e-13` absolute.

**This is reproduced deliberately and is not going to be fixed.** The angular
weights enter every moment and therefore every iteration count, and the entire
validation of this port rests on those counts matching the original bit for
bit. Correcting the quadrature would invalidate all of it.

The practical consequence: **`~1e-13` is the accuracy ceiling of every angular
integral in this solver.** Do not design a test that needs better.
`pybte.quadrature.leggauss_reference` gives the correct rule where a test
needs one.

## 5. GSIS cannot reach an arbitrarily tight tolerance in the ballistic limit

The iterate residual has a floor, `floor = beta x eps_macro`, where
`eps_macro` is the relative error of the HDG macroscopic solve. Measured at
`tau_R = 1e2` by varying `flow.tau_thr`, which scales `beta` directly:

| `tau_thr` | `beta_max` | residual floor | floor / beta |
|---|---|---|---|
| 1e+0 | 1.29e-03 | 3.15e-10 | 2.44e-07 |
| 1e-1 | 1.29e-04 | 3.11e-11 | 2.41e-07 |
| 1e-2 | 1.29e-05 | 3.16e-12 | 2.45e-07 |
| 1e-3 | 1.29e-06 | 3.17e-13 | 2.46e-07 |
| 1e-4 | 1.29e-07 | converges in 8 | — |

`eps_macro` grows as `Kn -> infinity` because the synthetic system degenerates
there: the momentum equation loses `q/Kn_R` and the stress equation loses
`N/Kn_C`, leaving a nearly singular operator. `cond_1(K)` goes `2.9e5` at
`tau_R = 1e-1`, `2.5e7` at `1`, **`2.5e11` at `1e2`**, where the macroscopic
temperature it returns is `|UQ_T| = 67` against a physical `|T| = 0.51`.
`beta` is what keeps that out of the answer.

**This is not a stall, and GSIS does not lose its convergence rate here.** At
any tolerance above the floor it takes *the same iteration count as CIS*,
which is correct when there is nothing to accelerate: 4/4 at `tol = 1e-6`,
5/5 at `1e-8` (`tau_R = 1e2`); 8/8 and 10/10 at `1e-8` and `1e-10`
(`tau_R = 1e1`).

**Practical rule: at `tau_R >= 10`, do not ask GSIS for better than `1e-8`.**
Loosen the tolerance, or run CIS, which needs only 4-13 iterations in this
regime. Raising `tau_thr` to force more acceleration makes it worse: with
`beta -> 1` the oversized macroscopic solution enters undamped and the answer
degrades from `3.3e-03` to `1.8e+00` as `tau_thr` goes `1 -> 1000`.

---

## Things that are not limitations, but surprise people

* **`Temp`, `Qx`, `Qy` are element integrals, not averages.** They carry the
  `INT_NODFUNC_TRI` weight, so `record.temp.sum()` is the domain integral of
  `T` — 0.25 for the shipped cavity, not the mean temperature. Divide by
  `mesh.tri_area` for cell averages.
* **`NAZIM` is silently forced even** (`nazim = (nazim // 2) * 2`).
* **The sweep is Gauss–Seidel, not Jacobi.** Within one direction, each
  element reads inflow from neighbours already updated in the same pass. This
  is what makes a sweep an exact solve for that direction, and it is why
  `TRI_ORDER`'s tie-breaking is reproduced exactly rather than replaced by any
  valid topological order.
* **CIS at `tau_R <= 1e-3` does not converge in any practical budget.** That
  is the phenomenon, not a defect. See #2 for why the residual will not tell
  you.
