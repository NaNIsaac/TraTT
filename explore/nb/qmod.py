"""Acoustic fields with a non-uniform local wavevector / phase (Entries 38+).

Phase field  phi(zeta, t) = S * A(s) * sin(theta(zeta, t)),  s = t - zeta (retarded time).
The uniform travelling wave of `Drive` is theta = q0 (zeta - t) = -q0 s. Generalisations:

  spatial_chirp(q0, beta, zc)   graded medium at fixed frequency: local wavevector
                                q(zeta) = q0 [1 + beta (zeta - zc)],  theta = int_0^zeta q - q0 t.
                                (The front is kept at zeta = t: the velocity grading is
                                neglected for the envelope. Approximation, stated.)
  retarded_chirp(q0, b)         frequency chirp at the transducer, non-dispersive medium:
                                theta = -Psi(s), Psi' = q0 (1 + b s): the local spatial
                                wavevector at depth zeta is q0 (1 + b (t - zeta)).
  q_wall(q1, q2, s_w, w)        retarded time s < s_w emitted at q1, later at q2 (continuous
                                phase); the wall sits at depth t - s_w, smoothed over w.
  phase_wall(q0, dphi, s_w, w)  same q0, carrier phase jumps by dphi for s > s_w.

`kinematic(alpha, d, eta, t)` is the first Born approximation of the same problem,
Dh = i pi alpha int_0^1 e^{-i phi(zeta)} e^{i eta (1 - zeta)} dzeta, evaluated with the full
e^{-i phi} (all harmonics), as the control for "is it just a Fourier transform".
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from tratt.drive import Drive


def _soft_step(x, w):
    return np.where(x >= 0, 1.0, 0.0) if w == 0 else 0.5 * (1.0 + np.tanh(x / w))


def _soft_ramp(x, w):
    """int of the soft step: (x)_+ for w = 0, smooth otherwise (continuous phase)."""
    if w == 0:
        return np.maximum(x, 0.0)
    return 0.5 * x + 0.5 * w * np.logaddexp(x / w, -x / w) - 0.5 * w * np.log(2.0)


def uniform(q0):
    return lambda z, t: q0 * (z - t)


def spatial_chirp(q0, beta, zc=0.5):
    return lambda z, t: q0 * (z + beta * (0.5 * z * z - zc * z)) - q0 * t


def retarded_chirp(q0, b):
    return lambda z, t: -q0 * ((t - z) + 0.5 * b * (t - z) ** 2)


def q_wall(q1, q2, s_w, w=0.0):
    return lambda z, t: -(q1 * (t - z) + (q2 - q1) * _soft_ramp((t - z) - s_w, w))


def phase_wall(q0, dphi, s_w, w=0.0):
    return lambda z, t: -q0 * (t - z) + dphi * _soft_step((t - z) - s_w, w)


@dataclass(frozen=True)
class PhaseDrive(Drive):
    """Drive with an arbitrary acoustic phase theta(zeta, t) and envelope A(s) (default step)."""

    theta_fn: Callable | None = None
    amp: Callable | None = None

    def A(self, s):
        s = np.asarray(s, dtype=np.float64)
        return self.amp(s) if self.amp else self.envelope(s)

    def phi(self, t):
        S = self.S
        th = self.theta_fn or uniform(self.q)

        def _phi(zeta):
            zeta = np.asarray(zeta, dtype=np.float64)
            return S * self.A(t - zeta) * np.sin(th(zeta, t))

        return _phi


def kinematic(alpha, d, eta, t, nz=8001):
    """First-Born |Dh|^2 on the eta grid at one time (vectorised, chunked over eta)."""
    z = np.linspace(0.0, 1.0, nz)
    w = np.full(nz, 1.0 / (nz - 1))
    w[[0, -1]] *= 0.5
    f = np.exp(-1j * d.phi(t)(z)) * w
    eta = np.asarray(eta, dtype=np.float64)
    out = np.empty(eta.size)
    for i in range(0, eta.size, 256):
        e = eta[i:i + 256, None]
        out[i:i + 256] = np.abs(np.pi * alpha * (np.exp(1j * e * (1 - z[None, :])) @ f)) ** 2
    return out


def features(x, y, rel=0.1):
    """Line-shape features of y(x) on a window: peaks above rel*max, FWHM of the envelope
    (outermost half-max crossings), centroid, skewness, central-minimum depth."""
    ymax = y.max()
    pk = [i for i in range(1, y.size - 1) if y[i] > y[i - 1] and y[i] >= y[i + 1] and y[i] > rel * ymax]
    above = np.flatnonzero(y > 0.5 * ymax)
    width = x[above[-1]] - x[above[0]] if above.size else np.nan
    p = np.trapezoid(y, x)
    c = np.trapezoid(x * y, x) / p if p > 0 else np.nan
    var = np.trapezoid((x - c) ** 2 * y, x) / p if p > 0 else np.nan
    skew = np.trapezoid((x - c) ** 3 * y, x) / p / var**1.5 if p > 0 else np.nan
    k0 = int(np.argmin(np.abs(x)))
    return dict(npk=len(pk), peaks=[float(x[i]) for i in pk], width=float(width), peak=float(ymax),
                integral=float(p), centroid=float(c), skew=float(skew), centre=float(y[k0]),
                xpeak=float(x[int(np.argmax(y))]))
