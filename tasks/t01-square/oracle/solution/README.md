# Phonon Boltzmann transport solver

Steady 2D2V gray linear Callaway phonon BTE.  Nodal discontinuous Galerkin in
space, discrete ordinates in angle, source iteration in time.

    from pybte import Case, Solver
    record = Solver(Case.from_yaml("cases/kn_1e-2.yaml")).run()
    record.iterations, record.converged, record.temp

or from the command line::

    python -m pybte cases/kn_1e-2.yaml

## The model

Phonons relax by two channels with relaxation times `tau_R` (resistive,
momentum-destroying) and `tau_N` (normal, momentum-conserving), combined as
`1/tau_C = 1/tau_R + 1/tau_N`.  With `Vg = L = 1` the Knudsen number is
`tau_R`: large means ballistic, small means diffusive, and in the diffusive
limit the temperature obeys Fourier's law with `kappa = Cv Vg^2 tau_R / 3`.
`pybte.analytic` carries that limit for the shipped cavity as a 200-term
series.

## How one iteration works

    while True:
        sweep()                 # exact transport solve, one direction at a time,
                                # in a topological order of the element graph
        compute_moments()       # T, qx, qy from the distribution
        if residual < tol: break

The sweep is Gauss-Seidel: within a direction each element reads inflow from
neighbours already updated in the same pass, which is what makes it an exact
solve for that direction.

## Cases

Two boundary-condition families on the same square and the same
discretisation (200 elements, `DEG = 2`, a `10 x 20` angular mesh,
`tol = 1e-8`), each at five `(tau_R, tau_N)` pairs spanning the diffusive,
transition, ballistic and hydrodynamic regimes:

    F1_*   all four walls isothermal, T = 1 on the north wall
    F2_*   east and west walls diffusely reflecting (adiabatic)

A diffusely reflecting wall re-emits, isotropically, the energy that hits
it, so its emission depends on the solution: `bc.py` recomputes it from the
distribution at the start of every iteration.

## Layout

    pybte/          the solver
      driver.py     the outer iteration
      sweep.py      the transport sweep
      moments.py    moments and the residual
      basis.py      nodal DG bases
      quadrature.py Gauss-Legendre and triangle rules
      bc.py         boundary conditions
      mesh/         gmsh reader, topology, geometry
      analytic.py   the diffusion-limit series
    cases/          case files
    meshes/         gmsh meshes
    tests/          pytest suite

## Notes

* `Temp`, `Qx`, `Qy` on the run record are element *integrals*, not averages:
  `record.temp.sum()` is the domain integral of `T`.  Divide by
  `mesh.tri_area` for cell averages.
* `NAZIM` is silently forced even.
* The angular quadrature is accurate to about `1e-13`, not to machine
  precision; do not design a test that needs better.
