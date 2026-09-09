"""Optional numba compilation.

``pybte`` runs correctly without numba -- every jitted kernel has a pure
Python/numpy twin and the test suite checks that the two agree.  When numba
*is* installed the sweep and the ordering pass use it, which is what makes
the per-iteration performance targets reachable.
"""
from __future__ import annotations

import os

__all__ = ["njit", "prange", "HAVE_NUMBA", "numba_available",
           "set_num_threads"]

_DISABLE = os.environ.get("PYBTE_NO_NUMBA", "").strip() not in ("", "0", "false", "False")

try:  # pragma: no cover - exercised implicitly
    if _DISABLE:
        raise ImportError("disabled by PYBTE_NO_NUMBA")
    import numba
    from numba import njit as _njit
    from numba import prange as _prange

    HAVE_NUMBA = True
    njit = _njit
    prange = _prange
except ImportError:  # pragma: no cover
    HAVE_NUMBA = False
    numba = None

    def njit(*args, **kwargs):
        """No-op stand-in: returns the undecorated function."""
        if len(args) == 1 and callable(args[0]) and not kwargs:
            return args[0]

        def deco(fn):
            return fn
        return deco

    prange = range


def numba_available() -> bool:
    return HAVE_NUMBA


def set_num_threads(n: int) -> int:
    """Cap the jitted kernels' thread count.  Returns what was actually set.

    ``n = 0`` means "leave numba's default alone".  Numba refuses to raise the
    count above ``NUMBA_NUM_THREADS``, so we clamp rather than let it raise.
    """
    if not HAVE_NUMBA or n <= 0:
        return 0
    n = min(int(n), numba.config.NUMBA_NUM_THREADS)
    numba.set_num_threads(n)
    return n
