# Krylov solution of the outer iteration

`scheme.method: krylov` replaces the *update rule* of source iteration, and
nothing else. The sweep, the discretisation, the boundary treatment and the
discrete fixed point are untouched.

## What it does

One outer step of source iteration is the affine map

```
u -> F(u) = M( A^-1 b(u) ) = T u + g
```

on the state `u = (T_s, qx_s, qy_s)` — build the equilibrium source from the
moments, invert the transport operator exactly by one sweep, take moments.
Repeating it is Richardson's method on `(I - T) u = g`, which converges at
the spectral radius of `T`; with energy-conserving collisions that is
`1 - O((l_mfp/L)^2)`, hence the stall documented in `LIMITATIONS.md` #2.

`krylov.py` hands the same system to GMRES: every matrix–vector product is
one sweep, the Krylov basis is grown without restarting (restarting throws
the space away and costs an order of magnitude here), orthogonalisation is
classical Gram–Schmidt applied twice, and the residual of every product is
known from the Givens-rotated Hessenberg matrix without an extra sweep.

**The fixed point is the same as source iteration's** — it is the solution
of the same linear system — reached to round-off rather than to a
tolerance. That is what makes this path usable as a reference generator,
which the shipped GSIS is not (`LIMITATIONS.md` #1).

## The state must hold everything the sweep reads back

With thermalising walls the state is the three moment fields alone: the
sweep order is topological, so one sweep is an exact block-triangular solve
and the result depends on the moments and the boundary inflow only. A
diffusely reflecting wall re-emits what it received, so its emission is
appended to the state; it is linear in the distribution and the map stays
affine. Two situations break exactness and make `run()` raise:

* a direction whose sweep graph had a cycle broken (part of the inflow is
  lagged);
* periodic faces (the partner element cannot be placed upwind; its outflow
  is lagged).

## Stopping and accounting

GMRES is pushed to `tol * krylov_tol_factor` (floored at round-off), with a
bounded extension once it is under `tol * 1e-2`. A closing sweep then turns
the converged state back into the distribution the run reports, and
**convergence is decided on the true transport residual** after that sweep
— never on the step between iterates. `iterations` counts sweeps: one for
`g`, one per Krylov product, one to close.

## Measured, shipped square, 200 elements, `DEG = 2`, `10 x 20`, `tol = 1e-8`

| `(Kn_R, Kn_N)` | source iteration | Krylov | transport residual |
|---|---|---|---|
| (0.001, 1e5) | > 200 000 | 992 | 2.2e-15 |
| (0.01, 1e5) | 16 830 | 228 | 2.8e-15 |
| (0.1, 1e5) | 317 | 37 | 8.2e-16 |
| (1, 1) | 33 | 21 | 3.0e-16 |
| (10, 0.01) | 2 706 | 574 | 7.5e-15 |

With diffusely reflecting east/west walls: 1 104 / 281 / 47 / 23 / 688
against > 200 000 / 32 343 / 643 / 75 / 2 444.

## Limits

* **Not Knudsen-independent.** The count still grows as the medium thickens
  (14 → 1 388 from `tau_R = 1` to `1e-5`), far more slowly than source
  iteration. The hydrodynamic corner `(10, 0.01)` is its weakest cell: 4–5x.
* **Memory.** The basis is `k` vectors of length `3 * NDOF * N_TRIS`
  (+ wall emission); `k ~ 1 400` at `tau_R = 1e-5` is 40 MB here and would
  be tens of GB on a production mesh. `krylov_max_bytes` caps it, at which
  point restarts begin and the count rises.
* **Preconditioning** by the macroscopic moment system (`krylov_precond`) is
  experimental and off: it helps in the diffusive limit and is harmful in
  the ballistic one, and no formulation valid across the whole range has
  been found.
* **Periodic faces** are not supported; see above.

## Configuration

```yaml
scheme:
  method: krylov          # default: source
  krylov_tol_factor: 1e-6
  krylov_max_bytes: 5e8
  krylov_max_rounds: 3
  krylov_precond: false
```
