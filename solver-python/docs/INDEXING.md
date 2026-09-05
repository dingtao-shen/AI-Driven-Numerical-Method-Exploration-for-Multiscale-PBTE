# Indexing conventions

The Fortran is 1-based and column-major; `pybte` is 0-based and (mostly)
row-major. This file is the single place the mapping is written down. Nothing
in the package mixes conventions: **every index in `pybte` is 0-based**, and
`-1` means "absent" wherever the Fortran uses `0`.

## Scalars

| Fortran | `pybte` |
|---|---|
| `NDOF_TRI` | `case.ndof_tri` = `(deg+1)(deg+2)/2` |
| `NDOF_FC` | `case.ndof_fc` = `deg+1` |
| `NP_TRI`, `NP_FC` | `case.np_tri`, `case.np_fc` |
| `N_TRIS`, `N_FCS` | `mesh.n_tris`, `mesh.n_faces` |
| `TAU_C` | `case.tau_c` = `1/(1/tau_r + 1/tau_n)` |

## The distribution function

The one place where the layout genuinely differs, because the hot loop wants
it to.

```
Fortran   VDF(L, TRID, J1, J2)        column-major
pybte     vdf[d, i, l]                C-contiguous,  d = j1 + npole*j2
```

Directions are flattened in Fortran order over `(J1, J2)`, so `d` runs fastest
over the polar index — the same order as `TRI_ORDER(:, J1, J2)` and
`cx.ravel(order="F")`. Putting `(direction, element)` first and the DOF last
makes the innermost solve contiguous; the Fortran layout would stride by
`N_TRIS` on every DOF access.

`Solver.vdf_fortran()` returns the Fortran-shaped view for stage comparison.

## Arrays kept in Fortran index order

These are compared elementwise against the dumps, so their axis order is
preserved exactly; only the base changes.

| Fortran | `pybte` | shape |
|---|---|---|
| `INT_NODFUNC_TRI(L,I)` | `integrals.int_tri[l, i]` | `(ndof_tri, n_tris)` |
| `INT_NODFUNC_TRI_TRI(M,L,I)` | `integrals.int_tri_tri[m, l, i]` | `(ndof_tri, ndof_tri, n_tris)` |
| `INT_NODFUNC_TRI_TRI_X(I,M,L)` | `integrals.int_tri_tri_x[i, m, l]` | `(n_tris, ndof_tri, ndof_tri)` |
| `INT_NODFUNC_TRI_FC(I,IL,M,L)` | `integrals.int_tri_fc[i, il, m, l]` | `(n_tris, 3, ndof_tri, ndof_tri)` |
| `INT_NODFUNC_TRI_FC_FC(M,L,I,IL)` | `integrals.int_tri_fc_fc[m, l, i, il]` | `(ndof_tri, ndof_fc, n_tris, 3)` |
| `INT_NODFUNC_FC_FC(M,L,F)` | `integrals.int_fc_fc[m, l, f]` | `(ndof_fc, ndof_fc, n_faces)` |
| `T_s(L,I)` | `mom.ts[l, i]` | `(ndof_tri, n_tris)` |
| `UQ(7*NDOF_TRI, I)` | `acc.uq[:, i]` | `(7*ndof_tri, n_tris)` |
| `U_TRACE(3*NDOF_FC, F)` | `acc.u_trace[:, f]` | `(3*ndof_fc, n_faces)` |

Note that `INT_NODFUNC_TRI_TRI_X` and `INT_NODFUNC_TRI_TRI` use *different*
axis orders in the Fortran itself — element index last in one, first in the
other. That is reproduced rather than tidied up, because the sweep indexes
both and the mass matrix is only symmetric to `1e-13`.

## Mesh topology

`Mesh` unpacks the Fortran's packed integer tables into named fields, and
`-1` replaces `0` for "no such entity".

| Fortran | `pybte` |
|---|---|
| `TRIANGLES_TAG(I,1:3)` node ids | `mesh.tri_nodes[i]`, 0-based |
| `TRIANGLES_TAG(I,4:6)` face ids | `mesh.tri_faces[i]`, 0-based |
| `TRIANGLES_INF(I,1)` area | `mesh.tri_area[i]` |
| `TRIANGLES_INF(I,2*IL:2*IL+1)` normal | `mesh.tri_normal[i, il]`, `il = IL-1` |
| `TRIANGLES_Hmin(I)` | `mesh.tri_hmin[i]` |
| `FACES_TAG(F,1:2)` nodes | `mesh.face_nodes[f]` |
| `FACES_TAG(F,3)` BC index, 0 = interior | `mesh.face_bc[f]`, **-1 = interior** |
| `FACES_TAG(F,4)` periodic pair / wall queue | `mesh.face_pair[f]`, -1 = none |
| `FACES_TAG(F,5:6)` triangle(+)/(-) | `mesh.face_tri[f, 0:2]`, -1 = none |
| `FACES_TAG(F,7:8)` local face id | `mesh.face_lfc[f, 0:2]`, 0-based |
| `FACES_INF(F)` length | `mesh.face_len[f]` |

`mesh.fortran_triangles_tag()`, `fortran_triangles_inf()` and
`fortran_faces_tag()` rebuild the packed 1-based tables for dump comparison.

### Local edge convention

Unchanged from `Spatial_Mesh.f90`:

```
   3
   | \        local edge 0 (Fortran 1):  node 0 -- node 1
   |  \       local edge 1 (Fortran 2):  node 1 -- node 2
   1---2      local edge 2 (Fortran 3):  node 2 -- node 0
```

`triangle(+)` is the element for which the face's stored node order runs
along the element's own edge direction; `triangle(-)` is the other one.

## The reference triangle and its bases

```
 (0,1)
   |\         element nodes:  (i/DEG, j/DEG) for j = 0..DEG, i = 0..DEG-j
   | \        stored with j outer, i inner
   |  \
   |___\      face nodes: -1 + 2k/DEG on the reference segment [-1, 1]
 (0,0) (1,0)
```

Monomials are grouped by total degree and, within a degree, by the `eta`
exponent, so `xi**a * eta**b` has index `(a+b)(a+b+1)/2 + b`
(`basis.monomial_index`). `nodfun_tri[m, k]` is the coefficient of monomial
`k` in shape function `m`.

## The GSIS block layout

`UQ` carries seven fields per element DOF and `U_TRACE` three per face DOF,
in this order:

```
UQ      block b, DOF m  ->  index b*ndof_tri + m
        b = 0 T,  1 qx,  2 qy,  3 Lxx,  4 Lxy,  5 Lyx,  6 Lyy

U_TRACE block b, DOF m  ->  index b*ndof_fc + m
        b = 0 T_hat,  1 qx_hat,  2 qy_hat
```

The global trace vector is `U_TRACE` flattened in Fortran order, i.e. face
index slowest: row `f*3*ndof_fc + b*ndof_fc + m`. This matches the Fortran's
`FFA(M + b*NDOF_FC + (I-1)*3*NDOF_FC)` exactly, so the assembled matrix can be
compared against the dumped CSR without permutation.
