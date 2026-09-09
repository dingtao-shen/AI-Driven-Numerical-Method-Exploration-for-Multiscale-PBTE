"""Angular (discrete-ordinate) mesh -- the port of ``Velocity_Mesh.f90``.

``THE`` is a Gauss-Legendre rule on ``[0, pi]`` **reversed**; ``PHI`` is two
Gauss-Legendre rules, on ``[0, pi]`` and ``[pi, 2pi]``, each reversed within
its own half.  The reversals are pure bookkeeping -- they do not change any
sum -- but they fix the direction ordering, and the direction ordering fixes
``TRI_ORDER``, so they are reproduced exactly.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constants import PI
from .quadrature import gauss_legendre

__all__ = ["VelocityMesh", "build_velocity_mesh"]


@dataclass(frozen=True)
class VelocityMesh:
    """Discrete ordinates and their solid-angle weights.

    Attributes
    ----------
    the, wthe : ndarray (npole,)
    phi, wphi : ndarray (nazim,)
    domega    : ndarray (npole, nazim)   ``sin(the)*wthe*wphi``; sums to 4*pi
    cx, cy    : ndarray (npole, nazim)   the in-plane group-velocity components
    """
    the: np.ndarray
    wthe: np.ndarray
    phi: np.ndarray
    wphi: np.ndarray
    domega: np.ndarray
    cx: np.ndarray
    cy: np.ndarray

    @property
    def npole(self) -> int:
        return self.the.size

    @property
    def nazim(self) -> int:
        return self.phi.size

    @property
    def ndir(self) -> int:
        return self.npole * self.nazim

    def flat(self):
        """(cx, cy, domega) flattened in Fortran order over (J1, J2)."""
        return (self.cx.ravel(order="F"),
                self.cy.ravel(order="F"),
                self.domega.ravel(order="F"))


def build_velocity_mesh(npole: int, nazim: int, vg: float = 1.0) -> VelocityMesh:
    nazim = (nazim // 2) * 2
    if nazim < 2:
        raise ValueError("nazim must be >= 2")
    half = nazim // 2

    a, w = gauss_legendre(npole, 0.0, PI)
    the = a[::-1].copy()
    wthe = w[::-1].copy()

    phi = np.zeros(nazim)
    wphi = np.zeros(nazim)

    a, w = gauss_legendre(half, 0.0, PI)
    phi[0:half] = a[::-1]
    wphi[0:half] = w[::-1]

    a, w = gauss_legendre(half, PI, 2.0 * PI)
    phi[half:nazim] = a[::-1]
    wphi[half:nazim] = w[::-1]

    sin_the = np.sin(the)
    domega = (sin_the * wthe)[:, None] * wphi[None, :]
    cx = np.repeat((vg * np.cos(the))[:, None], nazim, axis=1)
    cy = (vg * sin_the)[:, None] * np.cos(phi)[None, :]

    return VelocityMesh(the=the, wthe=wthe, phi=phi, wphi=wphi,
                        domega=domega, cx=cx, cy=cy)
