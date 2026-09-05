# Known limitations

What this solver cannot do, or does in a way you need to know about. Each
entry is a property of `pybte` **now** — none is a "to do".

Defects that were found in the Fortran original and *fixed* during the port
are not listed here; they no longer affect anything. They are recorded in the
git history at the `Port ACC_2D2V_LinearCallawayModel` commit and the ones
following it, together with the reference sources they were found in.

---

## 1. GSIS does not converge to the same answer as CIS

**The most important thing on this page.**

The two schemes are not two ways of solving the same discrete system. Both
converge, both to `residual_iterate < 1e-13`, and they land in different
places:

| | CIS | GSIS |
|---|---|---|
| iterations (`tau_R = 1e-1`) | 557 | 49 |
| transport residual at the answer | `4.3e-14` | `3.2e-03` |
| `int T dA` | 0.249930 | 0.249996 |

`||T_CIS - T_GSIS||_inf / ||T_CIS||_inf = 1.7e-2`.

The *transport residual* — `||A f - b||/||b||` for the discrete system,
evaluated at the converged state with no sweep — settles which one solves the
kinetic problem. CIS drives it to round-off, so CIS is at the discrete fixed
point. GSIS leaves it at `3e-3`, so GSIS is not.

The cause is structural. The macroscopic system is discretised by **HDG**, the
kinetic one by **upwind DG**, with different wall closures. GSIS is
fixed-point preserving when the macroscopic equations are the exact moment
equations of the kinetic discretisation; these are not. Take the converged CIS
solution, apply one macroscopic solve, and the recovered `UQ_T` differs from
the kinetic moment `T_VDF` by `1.7e-1` in the wall-adjacent elements and
`3e-3` in the bulk — the defect at its source.

**It does not go away under refinement.** `tools/fixed_point_study.py`:

| varying | values | gap |
|---|---|---|
| `DEG` = 1, 2, 3 | fixed `tau_R = 1e-1`, 200 elements | 1.60e-2, 1.65e-2, 1.69e-2 |
| elements = 50, 200, 800 | fixed `tau_R = 1e-1`, `DEG = 3` | 1.66e-2, 1.69e-2, 1.43e-2 |
| `tau_R` = 1, 1e-1, 1e-2 | fixed discretisation | 1.15e-2, 1.69e-2, 3.72e-3 |

Sixteen times the elements buys 14%. Only reducing `tau_R` closes the gap, as
the Knudsen layer thins and both schemes approach the same Fourier limit.

**What follows.** GSIS is a fast *approximate* solver, not an oracle for the
CIS answer. Anything that grades one scheme against the other will reject a
correct implementation. Grade instead against quantities that do not depend on
which scheme produced them:

* the transport residual (`iteration.true_residual: true`);
* the analytic Fourier limit at small `Kn` (`pybte.analytic`);
* the iteration counts, which are unaffected and are the real phenomenon —
  16 836 against 27 at `tau_R = 1e-2`.

Fixing it means rebuilding the macroscopic system as the exact moment system
of the kinetic discretisation. That is a redesign of the acceleration scheme,
not a repair, and it would change every GSIS iteration count.

## 2. Do not trust a CIS residual

The stopping criterion measures the **step** between iterates, not the
**error**. Source iteration converges linearly, `r_{n+1} = rho r_n`, so the
error remaining when you stop is the sum of all future steps,
`r rho/(1 - rho)` — and `rho -> 1` as the medium becomes optically thick.

| `tau_R` | CIS iterations | reported residual | actual error in `int T dA` |
|---|---|---|---|
| 1e-2 | 16 836 (converged) | 1.0e-08 | 1.0e-05 |
| 1e-3 | 200 000 (truncated) | 1.5e-06 | **4.4e-02** |
| 1e-4 | 200 000 (truncated) | 2.5e-06 | **2.2e-01**, i.e. 89% wrong |

Directly measured at `tau_R = 1e-2`, stopping at `tol = 1e-4`: 2 599
iterations, residual `9.99e-05`, actual error `1.46e-01` — **1465x larger** —
against a reference converged to `1e-12` in 31 245 iterations.

`RunRecord.error_estimate` fits `rho` to the tail of the residual history and
reports `r rho/(1-rho)`; it lands within an order of magnitude of the truth.
`RunRecord.contraction` gives `rho` itself. **Look at these before believing a
CIS result.**

Note that `iteration.true_residual: true` does *not* expose this. For CIS it
measures the same moment change the iterate residual does (`1.52e-06` against
`1.52e-06` on the truncated `tau_R = 1e-3` run). What it is for is #1.

The stopping rule is deliberately left as it was: it is the subject of the
benchmark this solver was ported to support.

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
