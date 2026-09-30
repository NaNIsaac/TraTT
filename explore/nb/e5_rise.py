"""Entry 5: finite transducer rise time.

With envelope A(s), s = t - zeta, depth zeta carries modulation S A(t - zeta).
If A varies slowly on the acoustic wavelength (rise >> 1/alpha_ac), harmonic n has
local coupling pi alpha J_n(S A), and at exact resonance eta = n q the Laue
rotations accumulate:

    I(n q, t) = sin^2( pi alpha Theta_n(t) ),  Theta_n(t) = | int_0^1 J_n(S A(t - zeta)) dzeta |

(for n = 0 the unmodulated part contributes J_0(0) = 1, which the integral
already includes). rise = 0 recovers Entry 1 and C-1.10.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

z = np.linspace(0.0, 1.0, 4001)
for S, alpha, rise in [(1.5, 2.0, 0.3), (1.5, 2.0, 0.8), (3.0, 2.0, 0.5), (1.5, 4.0, 0.5), (1.5, 2.0, 0.05)]:
    d = Drive(S=S, alpha_ac=20.0, rise=rise)
    q = d.q
    t = np.linspace(0.0, 1.0 + rise, 121) + 0.0037
    ns = (0, 1, 2)
    I = rocking_cube(alpha, d, np.array([n * q for n in ns]), t, n_z=1500)
    print(f"S={S} alpha={alpha} rise={rise}")
    for j, n in enumerate(ns):
        Th = np.array([abs(np.trapezoid(jv(n, S * d.envelope(tt - z)), z)) for tt in t])
        th = np.sin(np.pi * alpha * Th) ** 2
        # step-onset law for comparison: what you'd predict ignoring the rise
        Th0 = np.array([abs(np.trapezoid(jv(n, S * (tt - z >= 0)), z)) for tt in t])
        th0 = np.sin(np.pi * alpha * Th0) ** 2
        print(f"  n={n}: max|I-rise law|={np.max(np.abs(I[:, j]-th)):.3f}  "
              f"max|I-step law|={np.max(np.abs(I[:, j]-th0)):.3f}")
