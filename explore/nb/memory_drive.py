"""Arbitrary acoustic history for the transient-hunting entries (28+).

The crystal at slow time t holds the drive history s = t - zeta, zeta in [0, 1]: it is a
spatial record of the last t_tr of the transducer signal. Satellite n reads the
Bessel-weighted integral of that record,

    Theta_n(t) = int_0^1 J_n(S A(t - zeta)) * sign_n(t - zeta) dzeta,
    sign_n(s) = (-1)^n where the carrier phase is flipped by pi, else +1,

and the two geometries read Theta_n with different functions at exact resonance:
Laue sin^2(pi alpha |Theta_n|) (rotations about one axis add), Bragg
tanh^2(pi alpha |Theta_n|) (dR/dzeta = -i kappa (1 + R^2) gives R = -i tanh int kappa
for any real, depth-dependent kappa). Only phase jumps of 0 or pi are allowed: other
phases rotate about different axes and the angles no longer add.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.special import jv

from tratt.drive import Drive


@dataclass(frozen=True)
class HistoryDrive(Drive):
    """Drive with an arbitrary envelope A(s) and an optional pi carrier flip F(s) in {0, 1}."""

    env: Callable | None = None
    flip: Callable | None = None

    def envelope(self, s):
        s = np.asarray(s, dtype=np.float64)
        return np.where(s >= 0.0, self.env(np.clip(s, 0.0, None)), 0.0) if self.env else super().envelope(s)

    def phi(self, t):
        S, q = self.S, self.q

        def _phi(zeta):
            zeta = np.asarray(zeta, dtype=np.float64)
            s = t - zeta
            ph = np.pi * self.flip(s) if self.flip else 0.0
            return S * self.envelope(s) * np.sin(q * zeta - q * t + ph)

        return _phi


def theta(d: HistoryDrive, n: int, t, nz: int = 4001):
    """Signed accumulated coupling Theta_n(t) for each t."""
    z = np.linspace(0.0, 1.0, nz)
    out = []
    for x in np.atleast_1d(t):
        s = x - z
        w = jv(n, d.S * d.envelope(s))
        if d.flip:
            w = w * np.where(np.asarray(d.flip(s)) > 0.5, (-1.0) ** n, 1.0)
        out.append(np.trapezoid(w, z))
    return np.array(out)


def triangle(T):
    """0 -> 1 over T, 1 -> 0 over the next T, then 0."""
    return lambda s: np.clip(np.where(s < T, s / T, 2.0 - s / T), 0.0, 1.0)


def ramp(T):
    """Linear 0 -> 1 over T, then held."""
    return lambda s: np.clip(s / T, 0.0, 1.0)


def box(T):
    return lambda s: np.where(s < T, 1.0, 0.0)
