"""Non-ideal acoustic fields for the robustness study of the transient-only satellite (Entry 35).

Phase field, theta = q (zeta - t), s = t - zeta (retarded time):

    phi(zeta, t) = S * A(s, zeta) * sin(theta) + m2(s, zeta) * sin(2 theta + psi2)

A(s, zeta) may depend on depth (attenuation, dispersion) and m2 is a second-harmonic
phase index (minimal nonlinearity). With psi2 = 0 the Fourier coefficients of
e^{-i phi} at the satellite harmonic are real, the Laue rotations at resonance share
one axis, and the accumulated-coupling law

    Theta_n(t) = int_0^1 c_n(S A(t - zeta, zeta), m2(t - zeta, zeta)) dzeta,
    c_n(a, b) = (1/2pi) int e^{i (n th - a sin th - b sin 2th)} d th   (c_n(a, 0) = J_n(a))

holds; with psi2 != 0 it does not, and only full TT is used.

Envelopes (all of retarded time s, units of t_tr):
  ramp(T)                       linear 0 -> 1 over T (T = 0: step)
  dispersed_ramp(T, sig1)       ramp convolved with a Gaussian of width sig1 * sqrt(zeta)
  ringing_ramp(T, R, tau, f, psi)  ramp + R e^{-(s-T)/tau} sin(2 pi f (s - T) + psi) for s > T
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.special import erf

from tratt.drive import Drive

_TH = np.linspace(0.0, 2.0 * np.pi, 128, endpoint=False)


def c_n(n: int, a, b):
    """Fourier coefficient of e^{-i(a sin th + b sin 2th)} at e^{-i n th}; real for real a, b."""
    a = np.asarray(a, dtype=np.float64)[..., None]
    b = np.asarray(b, dtype=np.float64)[..., None]
    return np.mean(np.cos(n * _TH - a * np.sin(_TH) - b * np.sin(2 * _TH)), axis=-1)


def ramp(T: float) -> Callable:
    if T == 0.0:
        return lambda s, z=0.0: np.where(np.asarray(s) >= 0.0, 1.0, 0.0)
    return lambda s, z=0.0: np.clip(np.asarray(s, dtype=np.float64) / T, 0.0, 1.0)


def _ramp_conv(s, T, sig):
    """Linear ramp (0 -> 1 over T) convolved with a unit Gaussian of width sig (closed form)."""
    def F(x):  # (x)_+ convolved with the Gaussian
        u = x / sig
        return x * 0.5 * (1 + erf(u / np.sqrt(2))) + sig * np.exp(-0.5 * u * u) / np.sqrt(2 * np.pi)
    if T == 0.0:
        return 0.5 * (1 + erf(s / (sig * np.sqrt(2))))
    return (F(s) - F(s - T)) / T


def dispersed_ramp(T: float, sig1: float) -> Callable:
    """Envelope broadened by propagation: Gaussian width sig1 * sqrt(zeta) (diffusive spreading).
    Non-causal tails ahead of the front are kept (a dispersive precursor); A -> 1 behind."""
    def A(s, z):
        s = np.asarray(s, dtype=np.float64)
        sig = sig1 * np.sqrt(np.maximum(np.asarray(z, dtype=np.float64), 1e-12))
        sig = np.broadcast_to(sig, np.broadcast(s, sig).shape)
        tiny = sig < 1e-6
        out = np.where(tiny, ramp(T)(s), _ramp_conv(s, T, np.where(tiny, 1.0, sig)))
        return out
    return A


def ringing_ramp(T: float, R: float, tau: float, f: float, psi: float = 0.0) -> Callable:
    base = ramp(T)

    def A(s, z=0.0):
        s = np.asarray(s, dtype=np.float64)
        x = s - T
        ring = np.where(x > 0, R * np.exp(-np.clip(x, 0, None) / tau) * np.sin(2 * np.pi * f * x + psi), 0.0)
        return base(s) + ring
    return A


@dataclass(frozen=True)
class GeneralDrive(Drive):
    """Drive with depth-dependent envelope A(s, zeta), attenuation and a second harmonic.

    amp(s, zeta): envelope (default: the parent's step / smoothstep).
    gamma: amplitude attenuation per unit depth, A -> A e^{-gamma zeta}.
    eps2: second-harmonic index at the exit for full amplitude, m2 = eps2 zeta S (A e^{-gamma zeta})^2.
    psi2: phase of the second harmonic (0 keeps the law exact).
    """

    amp: Callable | None = None
    gamma: float = 0.0
    eps2: float = 0.0
    psi2: float = 0.0

    def A(self, s, z):
        s = np.asarray(s, dtype=np.float64)
        base = self.amp(s, z) if self.amp else self.envelope(s)  # envelopes handle s < 0 themselves
        return base * np.exp(-self.gamma * np.asarray(z, dtype=np.float64))

    def m2(self, s, z):
        return self.eps2 * np.asarray(z, dtype=np.float64) * self.S * self.A(s, z) ** 2

    def phi(self, t):
        S, q = self.S, self.q

        def _phi(zeta):
            zeta = np.asarray(zeta, dtype=np.float64)
            th = q * zeta - q * t
            s = t - zeta
            out = S * self.A(s, zeta) * np.sin(th)
            if self.eps2:
                out = out + self.m2(s, zeta) * np.sin(2 * th + self.psi2)
            return out

        return _phi


def theta_law(d: GeneralDrive, n: int, t, nz: int = 2001):
    """Signed accumulated coupling Theta_n(t) (psi2 = 0 only)."""
    z = np.linspace(0.0, 1.0, nz)
    out = []
    for x in np.atleast_1d(t):
        s = x - z
        # channel eta = +n q is driven by harmonic -n (cf. Entry 14, J_{-n}); with a second
        # harmonic c_{-n} != +-c_n, so the sign matters (checked against TT, Entry 35)
        out.append(np.trapezoid(c_n(-n, d.S * d.A(s, z), d.m2(s, z)), z))
    return np.array(out)
