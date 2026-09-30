"""Entry 15: coupling cancellation inside a ramp.

For S above the first zero of J_n (3.832 for n = 1), the ramp S A(s) sweeps J_n
through a sign change. The accumulated coupling Theta_n = |int J_n(S A) dzeta| then
rises and falls back while the ramp crosses the crystal, even if pi alpha Theta
stays < pi/2 -- a satellite that appears and vanishes with no Pendellosung.
Test with small alpha so sin^2 is monotone in Theta: I ~ sin^2(pi alpha Theta).
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

z = np.linspace(0.0, 1.0, 4001)
for S, alpha, rise, n in [(5.0, 0.5, 1.0, 1), (5.0, 0.5, 2.0, 1), (5.0, 0.5, 0.0, 1), (3.0, 0.5, 1.0, 1)]:
    d = Drive(S=S, alpha_ac=20.0, rise=rise)
    t = np.linspace(0.0, 1.0 + rise, 121) + 0.0037
    I = rocking_cube(alpha, d, np.array([n * d.q]), t, n_z=1500)[:, 0]
    Th = np.array([abs(np.trapezoid(jv(n, S * d.envelope(x - z)), z)) for x in t])
    th = np.sin(np.pi * alpha * Th) ** 2
    k = int(np.argmax(th))
    print(f"S={S} alpha={alpha} rise={rise} n={n}: max|I-law|={np.max(np.abs(I-th)):.3f} | "
          f"peak at t={t[k]:.2f} I={I[k]:.3f} (law {th[k]:.3f}) | final I={I[-1]:.3f} (law {th[-1]:.3f}) | "
          f"final/peak={I[-1]/I.max():.2f}")
