# Governing equations and discretisation

## 1. The physical model

Steady, two-dimensional, gray, **linearised Callaway** phonon Boltzmann
transport. The unknown is the perturbation `f(x, s)` of the phonon
distribution about a uniform reference, at position `x` and travel direction
`s` on the unit sphere.

Callaway's model splits phonon scattering into a **resistive** channel
(relaxation time `tau_R`, which destroys momentum and drives the system to a
local equilibrium at temperature `T`) and a **normal** channel (`tau_N`, which
conserves momentum and drives it to a *drifting* equilibrium). Writing
`tau_C = (1/tau_R + 1/tau_N)^-1` for the combined rate, the steady transport
equation is

```
c . grad f  +  f/tau_C  =  S
```

with the source collecting both channels,

```
S  =  Cv T_s / (4 pi tau_R)
   +  (1/tau_N) [ Cv T_s/(4 pi)  +  (3/(4 pi Vg^2)) (q_x c_x + q_y c_y) ]
```

Here `c = Vg s` is the group velocity, `Cv` the volumetric heat capacity and
`Vg` the group speed. The moments that close the system are

```
Cv T_s = <f>,        q_x = <c_x f>,        q_y = <c_y f>
```

with `<.>` the integral over the unit sphere. Because the model is linear,
`T_s` and `q` are *perturbations*: they may be negative, and the wall data are
perturbation temperatures.

The two limits worth keeping in view:

* `tau_R -> infinity` (ballistic): phonons stream from wall to wall.
* `tau_R -> 0` (diffusive): eliminating `f` gives Fourier's law with
  conductivity `kappa = Cv Vg^2 tau_R / 3`, so `T` satisfies Laplace's
  equation and `q = -kappa grad T`. With `Vg = 1` and a unit domain the
  Knudsen number is `Kn = tau_R`, which is why the shipped case at
  `tau_R = 1e-3` sits deep in the regime where source iteration crawls.
  `pybte.analytic` carries the Laplace solution for the shipped cavity.

### Boundary conditions

Incoming phonons at a wall (`c . n < 0`, `n` the outward normal) are set by
the wall:

* **thermalising** (`BC_TYP = 1`) — the wall emits an isotropic distribution
  at its own perturbation temperature, `f_w = Cv T_wall / (4 pi)`;
* **non-thermalising / adiabatic** (`BC_TYP = 2`) — diffusely reflecting: the
  wall emits isotropically at whatever level makes the net normal energy flux
  vanish,
  `f_w = -[integral over outgoing of (c.n) f] / [integral over incoming of |c.n|]`;
* **periodic** (`BC_TYP = 3`) — the incoming distribution is taken from the
  paired face, offset by `(BC_XOFF, BC_YOFF)`.

## 2. Angular discretisation — discrete ordinates

Directions are a tensor product of Gauss–Legendre rules in the polar angle
`theta` on `[0, pi]` and the azimuthal angle `phi` on `[0, pi]` and
`[pi, 2 pi]` (`NAZIM/2` nodes each, so `NAZIM` is forced even):

```
DOMEGA(j1, j2) = sin(theta_j1) w_theta(j1) w_phi(j2)      sums to 4 pi
c_x(j1, j2)    = Vg cos(theta_j1)
c_y(j1, j2)    = Vg sin(theta_j1) cos(phi_j2)
```

The problem is two-dimensional in space but the velocity space stays fully
three-dimensional — hence "2D2V": `c_z` never enters the transport operator,
but the polar integration over it is what produces the correct `1/3` in the
diffusion limit. Both `sum(DOMEGA c_x^2) = 4 pi Vg^2/3` and
`sum(DOMEGA) = 4 pi` are asserted in the unit tests.

## 3. Spatial discretisation — nodal DG

The domain is triangulated; on each element the solution is a degree-`DEG`
polynomial with `NDOF_TRI = (DEG+1)(DEG+2)/2` nodal degrees of freedom, and no
continuity is imposed between elements. Multiplying the transport equation by
a test function `phi_M` and integrating by parts over element `T` gives, for
one direction,

```
(1/tau_C) (phi_M, f)_T  -  (c . grad phi_M, f)_T  +  <phi_M, (c.n) f*>_dT  =  (phi_M, S)_T
```

with `f*` the **upwind** numerical flux: the element's own trace where
`c . n > 0`, the neighbour's (or the wall's) where `c . n < 0`. Splitting
`c.n = 0.5(c.n + |c.n|) + 0.5(c.n - |c.n|)` puts the outflow half on the left
and the inflow half on the right, which is exactly how `Solvers.f90` writes it:

```
A_SOL(M,L) = (1/tau_C) INT_NODFUNC_TRI_TRI(M,L,I)
           - c_x INT_NODFUNC_TRI_TRI_X(I,M,L) - c_y INT_NODFUNC_TRI_TRI_Y(I,M,L)
           + sum_IL 0.5 (speed + |speed|) INT_NODFUNC_TRI_FC(I,IL,M,L)
```

All the integrals are precomputed once (`pybte.integration`): the element ones
analytically, from monomial moments `int_T xi^a eta^b = a! b! / (a+b+2)!`
contracted with the nodal coefficients; the two that couple different
polynomial spaces or different elements by a 15-point Gauss–Legendre rule
along the face.

Because the inflow term is the *only* coupling between elements and it is
one-directional, the element systems can be solved in sequence rather than
simultaneously — that is the sweep.

## 4. The sweep

For each direction, elements are ordered so that every element's upwind
neighbours precede it (a topological sort of the per-direction dependency
graph, `pybte.ordering`). Visiting them in that order and solving each
`NDOF_TRI x NDOF_TRI` system in turn produces the exact solution of the whole
transport problem for that direction in one pass, because every inflow value
needed is already final.

This makes the sweep **Gauss–Seidel in the sweep ordering, not Jacobi**, and
the distinction is not cosmetic: a Jacobi port converges at a different rate,
so every iteration count downstream would drift. The ordering is reproduced
including its tie-breaking (ascending element index), so the port's iterate
sequence is bit-identical to the Fortran's.

`A_SOL` depends only on geometry, direction and `tau_C` — all iteration
invariant — so `pybte` LU-factorises it once at setup and reduces each sweep
to a triangular solve (128 MB on the shipped case at `DEG=3`).

## 5. CIS — conventional source iteration

```
repeat:
    f      <- sweep(T_s, q_s)          # transport, all directions
    T_s, q <- moments(f)               # angular integration
until  ||T - T_old|| / ||T||  <  TOL
```

The convergence rate is set by how much of the source is regenerated per
sweep. In the optically thick limit almost all of it is: the spectral radius
approaches 1 as `tau_R -> 0`, and the iteration count grows like `1/tau_R`.
Measured on the shipped mesh: 24 iterations at `tau_R = 1`, 317 at `1e-1`,
16836 at `1e-2`, and more than `2e4` (not converged) at `1e-3`. **This is the
phenomenon the whole project is about.**

## 6. GSIS — general synthetic iterative scheme

A macroscopic system is solved between sweeps and used to correct the
distribution. Its unknowns are seven fields per element DOF,
`UQ = [T, q_x, q_y, L_xx, L_xy, L_yx, L_yy]`, and three trace fields per face
DOF, `U_TRACE = [T_hat, q_hat_x, q_hat_y]`. The continuum content is

```
div q = 0                                   (steady energy balance)
q/tau_R + (Cv/3) grad T + div(HoT) = 0      (momentum, with a stress closure)
```

where the higher-order terms `HoT` are **computed from the kinetic solution**
rather than modelled — contracting `f` with the x- and y-derivative tensors
and the third-order angular moments `c_x(5c_x^2-3)`, `c_y(5c_x^2-1)`, and so
on, scaled by `tau_C/5`. That is what distinguishes a synthetic scheme from a
moment model.

The system is discretised by **HDG**: the element fields are eliminated
locally in favour of the face traces,

```
UQ_I = inv(AA_SOL_I) [ AA_SRC_I + sum_IL AA_TRACE_{I,IL} U_TRACE_face ]
```

leaving a sparse global system in the traces alone. That matrix depends only
on geometry, `tau_R`/`tau_C` and the stabilisation, so it is assembled and
factorised **once**; only the right-hand side changes per iteration. One
outer iteration is then:

```
f          <- sweep(T_s, q_s)
AA_SRC     <- HoT source from f
U_TRACE    <- global trace solve            (back-substitution only)
UQ         <- local recovery per element
f, T_s, q  <- correction
```

with the damped correction

```
tau_loc = tau_R / Hmin(I)
beta    = min(tau_loc, TAU_THR) / tau_loc
l_T     = (T_ACC - T_VDF) beta,     likewise l_qx, l_qy
f      += l_T Cv/(4 pi) + (c_x l_qx + c_y l_qy) (tau_C/tau_N) 3/(4 pi Vg^2)
```

`beta` damps the correction where the local element is optically thin, so the
acceleration cannot overshoot in the ballistic layer.

The payoff is that the iteration count stops depending on `tau_R`: 22, 30, 27,
50 iterations at `tau_R = 1, 1e-1, 1e-2, 1e-3` against CIS's 24, 317, 16836,
`>2e4`.

**But see `LIMITATIONS.md` #1**: because the macroscopic system is
discretised by HDG while the kinetic one is upwind DG, GSIS does *not*
converge to the same discrete fixed point as CIS. It is a fast approximate
solver, not an oracle for the CIS answer, and the two fields differ by
`~1.7e-2` at `tau_R = 1e-1`.

## 7. Convergence criterion

```
RESIDUAL = sqrt( sum_I (T_I - T_I^old)^2 / sum_I T_I^2 )
```

the normalised change in cell-average temperature between successive
iterates. Note what this is: a **step**, not an **error**.

Source iteration converges linearly, `r_{n+1} = rho r_n`, so the error left
after stopping is the sum of all future steps, `r rho/(1-rho)`. As the medium
becomes optically thick `rho -> 1` and that factor blows up. Measured at
`tau_R = 1e-2` on the shipped mesh, stopping at `tol = 1e-4`:

| reported residual | actual error against the converged answer |
|---|---|
| 9.99e-05 | 1.28e-01 relative in `int T dA`, against the same scheme at `tol = 1e-12` — **three orders of magnitude larger** |

That pseudo-convergence trap is preserved deliberately as the default
stopping rule, because it is one of the things the benchmark is meant to
probe. Two additional numbers are reported alongside it:

* `RunRecord.error_estimate` — `r rho/(1-rho)` with `rho` fitted to the tail
  of the residual history. This is the one that exposes the trap.
* `iteration.true_residual: true` — `||A f - b||/||b||` for the discrete
  transport system at the current state. For CIS this measures the same
  moment change the iterate residual does, so it does *not* expose the trap;
  what it does expose is that GSIS converges to a state which does not
  satisfy the transport system at all.
