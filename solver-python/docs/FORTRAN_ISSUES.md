# Defects and surprises in the Fortran reference

Everything here was found while porting `fortran-reference/ACC_2D2V_LinearCallawayModel`
and verified against the shipped case. Each entry says what the reference does,
what `pybte` does, and how to make the port reproduce the reference exactly if
you need to.

The reference tree is never edited. Where a fix is required to make the code
*run at all*, it is applied to a build copy by `tools/build_fortran.sh` and is
listed below.

| # | Severity | Where | Reproducible in `pybte`? |
|---|---|---|---|
| 1 | **critical** | `Synthetic_Acceleration.f90` | no — undefined behaviour |
| 2 | cosmetic | `Velocity_Distribution.f90` | `compat.qy_accumulation` |
| 3 | wrong output | `Out_Put_Result.f90` | `fortran_compat=True` |
| 4 | accuracy | `USD_Math.f90` | reproduced by default |
| 5 | **design** | `Synthetic_Acceleration.f90` | reproduced by default |
| 6 | **design** | `Synthetic_Acceleration1.f90` | `scheme.acc_variant: B` |
| 7 | correctness | `Solvers.f90` | see §7.1 of the proposal |
| 8 | correctness | `Synthetic_Acceleration.f90` | fixed; no switch |
| 9 | **correctness** | `Synthetic_Acceleration.f90` | fixed; no switch |
| 10 | determinism | `Velocity_Distribution.f90` | `restart.error_on_stale` |
| 11 | performance | `Synthetic_Acceleration.f90` | fixed; no switch |
| 12 | robustness | `Matrix.f90` | `scheme.on_cycle` |

---

## 1. `AA_TMP` is used uninitialised — the GSIS path is meaningless without a fix

`Init_Acceleration_New` allocates the face block

```fortran
ALLOCATE(AA(NDOF_FC,NDOF_FC), AA_TMP(3*NDOF_FC,3*NDOF_FC))
```

and then only ever assigns its three diagonal `NDOF_FC`-blocks:

```fortran
AA_TMP(        1:  NDOF_FC,        1:  NDOF_FC) = AA
AA_TMP(NDOF_FC+1:2*NDOF_FC,NDOF_FC+1:2*NDOF_FC) = AA
AA_TMP(2*NDOF_FC+1:3*NDOF_FC,2*NDOF_FC+1:3*NDOF_FC) = AA
```

The six off-diagonal blocks are never written. Fortran leaves them holding
whatever the allocator returned — in practice the just-freed `AA_SOL`
workspace, entries of order `1e5`. `AA_TMP` is copied into the diagonal block
of **every row** of the HDG global matrix, so the entire trace system is
corrupted.

Observed with gfortran 13.3 on the shipped case: the GSIS run "converges" in
5 iterations to a temperature field of magnitude `1e-27`, against a correct
mean of `0.25`. The matrix row for face 10 contains entries of `1e5` next to
legitimate entries of `1e-4`.

This is undefined behaviour, not a numerical choice, so there is nothing to
port faithfully. `tools/build_fortran.sh` inserts `AA_TMP = 0.d0` into its
build copy (`FIX-1`); `--no-fixes` reproduces the broken build for comparison.
With the fix, GSIS converges in 50 iterations to `mass = 0.2499975`.

## 2. `Qy` is never zeroed between iterations

`Calculate_Macro_Properties` opens with

```fortran
Temp = 0.d0
Qx = 0.d0
Qx = 0.d0     ! <- Qy intended
```

so the cell-average `Qy` accumulates across every iteration. Nothing in the
Fortran reads `Qy` — it is written to `RunTime.txt` by nobody and to the field
file from `VDF` directly — so the defect is invisible there. It is visible in
the stage dumps, which is how it was found.

`pybte` zeroes it. `compat.qy_accumulation: true` reproduces the reference
array exactly, for stage comparison.

## 3. The analytic reference heat flux is a factor of two too small

`Out_Put_Result.f90` builds the Laplace series

```
T = (2/pi) * sum_m [((-1)^(m+1)+1)/m] sin(m pi x) sinh(m pi y)/sinh(m pi)
```

then differentiates it to get the flux — but applies only `-Cv/3*TAU_R` to
the raw sum, dropping the `2/pi * m*pi = 2m` that survives differentiation.
The `Conduction_A.dat` fluxes are therefore exactly half the correct
`q = -(Cv Vg^2 tau_R/3) grad T`. The temperature column is correct.

`pybte.analytic.fourier_flux` is correct by default;
`fortran_compat=True` reproduces `Conduction_A.dat`.

## 4. `GaussLegendre` converges linearly, not quadratically

`USD_Math.f90` forms the Legendre derivative as

```fortran
Lp = N2*( L(:,N1) - y*L(:,N2) ) / (1.0d0-y*y)
```

with `N2 = n+1` where the identity `P_n'(y) = n(P_{n-1} - y P_n)/(1-y^2)`
calls for `n`. The weight formula carries the same `N2**2/N1**2` factor, so
the weights stay *consistent with whatever roots come out* — but the Newton
step is scaled by `n/(n+1)`, which degrades quadratic convergence to linear
with ratio `1/(n+1)`. Combined with a step-size exit test (`|dy| < 1e-13`)
the roots stop about `1e-13/n` short of the true zeros.

Measured: nodes differ from `numpy.polynomial.legendre.leggauss` by up to
`1.8e-14` (n=2) to `4.4e-16` (n=20); `sum(w)` on `[-1,1]` is short of 2 by
`~2e-13`. Downstream, `sum(DOMEGA)` misses `4*pi` by `3e-13` absolute.

`pybte` reproduces the iteration verbatim — the angular weights feed every
moment and hence every iteration count, so "fixing" it would break the
comparison. `pybte.quadrature.leggauss_reference` gives the correct rule for
tests.

## 5. GSIS does not preserve the kinetic fixed point

**This is the most consequential finding, and it contradicts the assumption in
Proposal 1 §10 that "CIS and GSIS converge to the same discrete fixed point to
`rtol=1e-8` — the load-bearing property for the whole benchmark".**

Measured on the shipped mesh at `DEG=3`, both schemes converged to
`residual_iterate < 1e-13`:

| | CIS | GSIS |
|---|---|---|
| iterations (`tau_R=1e-1`) | 557 | 49 |
| true transport residual | `4.3e-14` | `3.2e-03` |
| `int T dA` | 0.249930 | 0.249996 |

`||T_CIS - T_GSIS||_inf / ||T_CIS||_inf = 1.69e-2`, in the Fortran as well as
in the port (Fortran CIS 317 iterations, GSIS 30, gap `1.687e-2`; `pybte`
reproduces both — CIS bit-identically, GSIS to `5e-15`).

The *true transport residual* — `||A f - b||/||b||` of the discrete system,
evaluated with no sweep — settles it: CIS drives it to round-off, so CIS is at
the discrete kinetic fixed point. GSIS leaves it at `3e-3`, so GSIS is not.

The mechanism is visible in one step. Take the converged CIS solution, apply a
single GSIS macroscopic solve, and compare the recovered `UQ_T` against the
kinetic moment `T_VDF`. If GSIS were fixed-point preserving these would be
equal; instead they differ by `1.7e-1` in the wall-adjacent elements and
`3e-3` in the bulk (`tools/fixed_point_study.py` reports both).

The cause is structural, not a coding slip: the macroscopic system is
discretised by **HDG**, the kinetic system by **upwind DG**. They are two
different discretisations of the same continuum equations, with different wall
closures, so their solutions differ — most in the Knudsen layer, where the two
wall treatments disagree most. GSIS as published is fixed-point preserving
when the macroscopic equations reuse the kinetic discretisation; this
implementation does not.

`tools/fixed_point_study.py` measures how the gap scales. `gap` is
`||T_CIS - T_GSIS||_inf/||T_CIS||_inf`, `consist` the one-step defect above,
`bulk` its median over elements more than 0.2 from any wall:

| `tau_R` | `DEG` | `N_TRIS` | CIS it | CIS true res | GSIS it | GSIS true res | gap | consist | bulk |
|---|---|---|---|---|---|---|---|---|---|
| 1e+00 | 3 | 200 | 38 | 2.5e-15 | 35 | 6.0e-03 | 1.15e-02 | 1.13e+00 | 6.4e-02 |
| 1e-01 | 3 | 200 | 557 | 4.3e-14 | 49 | 3.2e-03 | 1.69e-02 | 1.66e-01 | 3.2e-03 |
| 1e-02 | 3 | 200 | 34846 | 9.6e-14 | 43 | 2.3e-04 | 3.72e-03 | 4.1e-02 | 5.0e-04 |
| 1e-01 | 1 | 200 | 556 | 3.4e-14 | 50 | 3.5e-03 | 1.60e-02 | 7.6e-02 | 4.5e-03 |
| 1e-01 | 2 | 200 | 557 | 3.8e-14 | 49 | 2.7e-03 | 1.65e-02 | 8.4e-02 | 3.2e-03 |
| 1e-01 | 3 | 50 | 557 | 6.5e-14 | 44 | 4.0e-03 | 1.66e-02 | 6.2e-02 | 3.5e-03 |
| 1e-01 | 3 | 800 | 557 | 2.4e-14 | 74 | 1.3e-03 | 1.43e-02 | 6.1e-01 | 2.8e-03 |

Read the three blocks separately.

* **Knudsen number** (rows 1–3, fixed discretisation): the gap falls with
  `tau_R` — 1.2e-2, 1.7e-2, 3.7e-3 — as the Knudsen layer thins and both
  schemes approach the same Fourier limit.
* **Polynomial order** (rows 2, 4, 5, fixed `tau_R` and mesh): 1.60e-2,
  1.65e-2, 1.69e-2 for `DEG` = 1, 2, 3. The gap does **not** fall; it drifts
  very slightly upward.
* **Mesh refinement** (rows 2, 6, 7, fixed `tau_R` and order): 1.66e-2,
  1.69e-2, 1.43e-2 for 50, 200 and 800 elements — a 16-fold increase in
  element count buys about 14%.

So this is *not* a discretisation error that vanishes under refinement. At
fixed physics the two schemes converge to two genuinely different answers, and
refining the grid does not bring them together; only reducing `tau_R` does.

Note also that CIS's iteration count is set by optical thickness alone: 557
iterations at `tau_R = 1e-1` regardless of `DEG` (1, 2, 3) or mesh (50, 200,
800 elements).

Consequences for a benchmark built on this solver:

* GSIS is a **fast approximate** solver, not an oracle for the CIS answer.
  A task graded on "does GSIS reproduce the CIS field" fails by `1.7e-2` on
  a correct implementation.
* The defensible oracle is the *true transport residual*, which is scheme
  independent, together with the analytic Fourier limit at small `Kn`.
* The iteration-count contrast (CIS 557 vs GSIS 49, growing to 16836 vs ~50 at
  `tau_R=1e-2`) is real and is unaffected by this.

## 6. Acceleration variant B diverges at every tested Knudsen number

`Synthetic_Acceleration1.f90` rescales the momentum and stress equations by
`TAU_R`, which is algebraically the same system as variant A but far better
conditioned as `TAU_R -> 0`. It is not the variant the Makefile builds.

Built and run (`tools/build_fortran.sh --variant B`), it diverges
geometrically from the first iteration at `tau_R = 1e-3`:
`mass` goes `5.8 -> -193 -> 6772 -> -2.4e5 -> ...`, reaching `NaN` by
iteration ~30 and then spinning to `TMAX`. `pybte` with
`scheme.acc_variant: B` reproduces this, and the divergence is not a
conditioning artefact — the 1-norm condition estimates of the two global
matrices are comparable:

| `tau_R` | cond₁(K) variant A | cond₁(K) variant B | A iters | B |
|---|---|---|---|---|
| 1e+00 | 2.5e7 | 4.9e7 | 22 | diverges |
| 1e-01 | 2.8e5 | 6.2e4 | 30 | diverges |
| 1e-02 | 1.6e4 | 1.4e5 | 27 | diverges |
| 1e-03 | 3.2e4 | 4.0e6 | 50 | diverges |

Variant B also drops the thermalising-wall rows from `BA_SOL` entirely (they
are commented out), so all three trace fields must be imposed from the DVM at
a wall; `pybte` does that, which is why its stage arrays match the variant-B
Fortran to `~3e-12`. The divergence is reproduced faithfully, not worked
around: it is a property of the shipped variant.

## 7. The boundary-condition dispatch is commented out

In `Solvers.f90` the whole `BC_TYP` branch is disabled and the *thermalising*
formula is applied unconditionally to every boundary face:

```fortran
!IF (BC_TYP(BCID).EQ.1) THEN
     FW = Cv/4.d0/PI*BC_TEMP(BCID)*INT_NODFUNC_TRI_TRI_FC(TRID,IL,:,1)
     A_SRC = A_SRC - 0.5d0*(speed-ABS(speed))*FW
!END IF
```

The diffusely-reflecting machinery exists in `Calculate_FLUX_WALL` but is
never called; the periodic branch is dead even though `Spatial_Mesh.f90` does
the face pairing. Every shipped case uses type 1 only, so the reference is
self-consistent — but `BC_TYP = 2` or `3` in a `control.in` silently produces
a hot wall.

`pybte` implements all three and dispatches on the type (proposal §7.1).
`BC_TYP = 4` (symmetry) is carried through the mesh reader but never
implemented anywhere; `pybte` raises rather than silently treating it as a
wall.

## 8. Periodic HDG assembly ignores the mirrored trace basis

Found while implementing §7.1, and latent in the reference because its
periodic branch never runs.

A face's 1-D trace basis is parameterised from its own node 1 to its own node
2. A periodic pair is traversed in *opposite* senses — both boundaries are
walked the same way around the domain — so DOF `k` of one face and DOF `k` of
its partner sit at mirrored points. The HDG assembly must permute
`k -> NDOF_FC-1-k` on contributions crossing the pair. Without it, GSIS on a
periodic channel converges to a field with a spurious `5e-2` variation along
the periodic direction and `|qx| ~ 3e-3` where symmetry demands zero, while
CIS on the same mesh gives `T` constant along `x` to `4e-8`.

`pybte` detects the orientation geometrically (`bc.periodic_flip`) and applies
the permutation in both the matrix and the right-hand side. After the fix the
`x`-variation drops to `1.1e-5` at `tau_R = 1e-2`.

## 9. The adiabatic wall's macroscopic trace flux is not projected

Found while testing §7.1, and latent in the reference for the same reason as
#8: `Solvers.f90` treats every wall as thermalising, so `BC_TYP = 2` never
occurs in a working run and this path is never executed.

An adiabatic wall is defined by `n . q = 0`. `Global_Problem_Solver_ACC`
imposes the macroscopic trace heat flux there from the kinetic solution, and
the tangential projection that would enforce the condition is written out --
but commented out in variant A:

```fortran
FFA(M+  NDOF_FC+(I-1)*3*NDOF_FC) = qxbc !qxbc*ny*ny - qybc*nx*ny
FFA(M+2*NDOF_FC+(I-1)*3*NDOF_FC) = qybc !-qxbc*nx*ny + qybc*nx*nx
```

Variant B has the same two lines with the projection *live*. So variant A
imposes the raw DVM flux, a trace with a non-zero normal component at a wall
that is by definition adiabatic.

The inconsistency is not cosmetic. On the adiabatic cavity at `tau_R = 1e-1`:

| | iterations | converged | `int T dA` |
|---|---|---|---|
| raw flux, as variant A ships | 300 (capped) | no | -3.8e+77 |
| tangential projection restored | 30 | yes | 0.5000000 |

0.5 is the exact answer: with `T = 0` on one wall, `T = 1` on the opposite one
and the sides adiabatic, the problem is one-dimensional and antisymmetric.

`pybte` applies the projection `q_hat = (I - n n^T) q` unconditionally, for
both variants. Nothing is lost by doing so, because there is no reference
behaviour to reproduce -- the branch is unreachable in the Fortran. The
kinetic side needs no such fix: CIS with adiabatic walls converges on its own
to `int T dA = 0.49995`, and the net normal flux through each adiabatic face
vanishes to round-off.

## 10. The restart file is read silently

`Init_Velocity_Distribution_Function` opens `VDF<P..T..NP..NA..>.out` and, if
it exists, uses it as the initial condition with only a line on stdout. A
stale file left in the working directory changes the initial VDF and therefore
the iteration count, with no other warning — a determinism hazard for anything
that grades on iteration counts.

`pybte` requires restart to be configured explicitly, records the file's
SHA-256 in the run record, and **errors** if a restart file is present while
restart is disabled (proposal §7.3).

## 11. Global assembly is quadratic in face count

`Init_Acceleration_New` zeroes the dense row-block workspace
`KKA(3*NDOF_FC, N_FCS*3*NDOF_FC)` inside the loop over faces, and then scans
all `N_FCS*3*NDOF_FC` columns of each row to compress it. Memory is fine, but
both the zeroing and the scan make assembly `O(N_FCS^2)`.

`pybte` accumulates each face's row block into a small dict of column blocks
and emits COO triplets directly, which is linear (proposal §7.4). The
resulting matrix matches the Fortran CSR to `1.5e-15` relative with identical
sparsity (213 760 nonzeros on the shipped case).

## 12. `Init_Triangle_Order` has no cycle detection

The topological sort loops

```fortran
DO
   DO I = 1, N_TRIS
      ... emit any element with no remaining incoming face ...
   END DO
   IF (NUM.EQ.N_TRIS) EXIT
END DO
```

with no test for a pass that emits nothing. On a mesh where some direction
produces a cyclic element dependency the loop spins forever, silently.

`pybte` raises `SweepCycleError` by default and offers
`scheme.on_cycle: break`, which emits the element with the fewest incoming
faces (lagging its inflow to the previous iterate) and records how many
directions needed it in the run record (proposal §7.7).

---

## Not defects, but worth knowing

* **`Temp`, `Qx`, `Qy` are element integrals, not averages.** They are
  weighted by `INT_NODFUNC_TRI`, so `SUM(Temp)` is the domain integral of `T`
  — `0.25` for the shipped cavity, not the mean temperature. `MASS` in the
  iteration log is that integral.
* **The residual is the change between iterates** -- a *step*, not an
  *error*. For a linearly converging sequence the remaining error is the sum
  of all future steps, `r * rho/(1-rho)`, and `rho -> 1` as the medium becomes
  optically thick. The trap this sets is severe and is preserved deliberately:

  | `tau_R` | CIS iterations | reported residual | actual error in `int T dA` |
  |---|---|---|---|
  | 1e-2 | 16 836 (converged) | 1.0e-08 | 1.0e-05 |
  | 1e-3 | 200 000 (truncated) | 1.5e-06 | **4.4e-02** |
  | 1e-4 | 200 000 (truncated) | 2.5e-06 | **2.2e-01**, i.e. 89% wrong |

  A residual of `2.5e-6` after 200 000 iterations looks like a nearly
  converged run. It is not. `RunRecord.error_estimate` fits `rho` to the tail
  of the residual history and reports `r rho/(1-rho)`, which lands within an
  order of magnitude of the actual error.

  Note that `iteration.true_residual: true` does **not** expose this. It
  reports the transport residual of the discrete system, which for CIS
  measures the same moment change the iterate residual does (1.52e-06 against
  1.52e-06 on the truncated `tau_R = 1e-3` run). What it *is* good for is
  separating the schemes: see #5.
* **`NAZIM` is silently forced even** (`NAZIM = (NAZIM/2)*2`).
* **The sweep is Gauss–Seidel, not Jacobi.** `VDF(:,TRIDext,...)` is read
  after upwind neighbours have already been overwritten in the same pass. A
  Jacobi port converges at a different rate and would invalidate every
  iteration count.
* **`A_SOL` is rebuilt and re-factorised every element, every direction, every
  iteration** although it depends only on geometry, direction and `TAU_C`.
  `pybte` hoists it (128 MB on the shipped case) and falls back to on-the-fly
  factorisation via `performance.precompute_inverse: false`.
