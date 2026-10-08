"""Entry 63: how a surface acoustic wave relates to the filled-depth picture.

A SAW is NOT this model's geometry: its wavevector lies in the surface and its
amplitude decays into the depth, so the satellites are displaced in an in-plane
direction and a one-dimensional depth solver cannot represent it. What the
solver CAN isolate is the half of the difference that is about depth: replace
the filled-depth front by a static exponential skin of depth d and ask what the
accumulated-coupling law reads.

    front (this work)   A(zeta) = 1 for zeta < t          grows with time
    skin  (SAW-like)    A(zeta) = exp(-zeta L / d)        fixed by frequency

Run: python explore/nb/e63_saw.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from tratt.core import rocking  # noqa: E402
from tratt.drive import Drive  # noqa: E402


def theta(n: int, S: float, zeta_d: float, n_pts: int = 20001) -> float:
    """Accumulated coupling of a static skin, Theta_n = int_0^1 J_n(S A(zeta)) d zeta."""
    z = np.linspace(0.0, 1.0, n_pts)
    return float(np.trapezoid(jv(n, S * np.exp(-z / zeta_d)), z))


def skin_phi(S: float, q: float, zeta_d: float, t: float):
    """phi(zeta) for a static exponential skin carried by a travelling phase."""

    def _phi(zeta):
        zeta = np.asarray(zeta, dtype=np.float64)
        return S * np.exp(-zeta / zeta_d) * np.sin(q * zeta - q * t)

    return _phi


def law_vs_tt():
    """Does the accumulated-coupling law still hold when the drive is a surface skin?"""
    alpha, alpha_ac, n = 3.0, 20.0, 1
    d = Drive(S=1.0, alpha_ac=alpha_ac)
    n_z = 4096
    print(f"Laue, alpha={alpha}, alpha_ac={alpha_ac}, n={n}, q={d.q:.2f}, n_z={n_z}")
    print(f"{'S':>5} {'zeta_d':>7} {'Theta_n':>9} {'law sin^2':>10} {'full TT':>9} {'diff':>9}")
    worst = 0.0
    for S in (0.5, 1.5, 3.0):
        for zeta_d in (0.02, 0.05, 0.1, 0.3, 1.0):
            th = theta(n, S, zeta_d)
            law = np.sin(np.pi * alpha * abs(th)) ** 2
            # average over acoustic phase: the skin is static, the carrier still moves
            ts = np.linspace(0.0, d.period, 8, endpoint=False)
            tt = float(np.mean([rocking(alpha, d.eta_channel(n), phi=skin_phi(S, d.q, zeta_d, t),
                                        n_z=n_z) for t in ts]))
            worst = max(worst, abs(tt - law))
            print(f"{S:5.1f} {zeta_d:7.3f} {th:9.5f} {law:10.5f} {tt:9.5f} {tt - law:9.2e}")
    print(f"worst |TT - law| = {worst:.2e}")


def effective_depth():
    """A skin of depth d reads like a filled depth of what?  Theta_n(skin) = Theta_n(front)."""
    print("\neffective filled depth of an exponential skin (n = 1)")
    print(f"{'S':>5} {'zeta_d':>7} {'Theta_skin':>11} {'t_eff':>7} {'t_eff/zeta_d':>13}")
    for S in (0.5, 1.5, 3.0):
        for zeta_d in (0.02, 0.05, 0.1, 0.3):
            th = theta(1, S, zeta_d)
            # a filled front of depth t gives Theta = t * J_1(S)
            t_eff = th / jv(1, S)
            print(f"{S:5.1f} {zeta_d:7.3f} {th:11.6f} {t_eff:7.4f} {t_eff / zeta_d:13.3f}")


def saw_bound():
    """The frequency- and material-independent cap on how much of Lambda_P a SAW can modulate."""
    rho, c11, c44 = 2329.0, 165.7e9, 79.6e9
    cL, cT = np.sqrt(c11 / rho), np.sqrt(c44 / rho)
    print()
    print("SAW depth bound (Si).")
    print(f"  c_L = {cL:.0f} m/s, c_T = {cT:.0f} m/s")
    for cR in (4920.0, 5378.0):
        kT = np.sqrt(1 - (cR / cT) ** 2)          # slower-decaying Rayleigh partial wave
        dec = 1.0 / (2 * np.pi * kT)              # decay length in units of Lambda_ac
        print(f"  c_R = {cR:.0f} m/s -> kappa_T/q = {kT:.3f}, skin = {dec:.3f} Lambda_ac, "
              f"modulated fraction of Lambda_P = {dec:.3f} R")
    dec = 1.0 / (2 * np.pi * np.sqrt(1 - (4920.0 / cT) ** 2))
    print()
    print(f"  resolved channels need R <= 0.34-0.57 (Entry 54), so the fraction is capped at "
          f"{dec * 0.34:.3f}-{dec * 0.57:.3f}")

    print()
    print(f"{'f GHz':>7} {'Lam_ac um':>10} {'skin um':>8} | "
          f"{'Si440 17keV R':>14} {'frac':>6} | {'Si440 30keV R':>14} {'frac':>6}")
    for f in (0.2, 0.5, 1.0, 2.0, 5.0):
        lam = 4920.0 / (f * 1e9) * 1e6
        skin = dec * lam
        out = f"{f:7.2f} {lam:10.3f} {skin:8.3f} |"
        for LamP in (54.34, 98.76):
            out += f" {lam / LamP:14.4f} {skin / LamP:6.3f} |"
        print(out)


if __name__ == "__main__":
    law_vs_tt()
    effective_depth()
    saw_bound()
