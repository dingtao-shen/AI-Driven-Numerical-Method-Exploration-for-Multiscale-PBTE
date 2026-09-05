"""Element recovery -- ``Local_Problem_Solver_ACC``.

    UQ(:, I) = AA_SRC(:, I) + sum_{IL=1..3} AA_TRACE(:, :, IL, I) @ U_TRACE(:, face(I, IL))

``AA_SRC`` and ``AA_TRACE`` are both already premultiplied by ``inv(AA_SOL)``
during setup, so this is three small matrix-vector products per element and
nothing else.
"""
from __future__ import annotations

import numpy as np

__all__ = ["local_solve"]


def local_solve(aa_src, aa_trace, u_trace, tri_faces, out=None):
    """``aa_src`` (7nd, n_tris); ``aa_trace`` (n_tris, 3, 7nd, 3nf);
    ``u_trace`` (3nf, n_faces); returns ``UQ`` (7nd, n_tris)."""
    if out is None:
        uq = aa_src.copy()
    else:
        uq = out
        uq[:] = aa_src
    traces = u_trace[:, tri_faces]                       # (3nf, n_tris, 3)
    uq += np.einsum("ilpr,ril->pi", aa_trace, traces, optimize=True)
    return uq
