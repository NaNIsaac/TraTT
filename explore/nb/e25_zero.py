"""Entry 25: approaching the Bessel zero of J_1 (S -> 3.8317).

Satellite-only law: I(q, t) = sin^2(pi alpha J_1(S) t) -> 0 and tau_1 -> infinity.
What survives at the zero comes from couplings the law drops (other harmonics,
main-line tail). Measure: (a) centre I at t = 1 vs S across the zero, (b) the
residual at the zero vs alpha_ac and alpha, (c) where the minimum sits in S.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros, jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

Z = jn_zeros(1, 1)[0]
t = np.linspace(0.0, 1.0, 51) + 0.0037


def centre(S, alpha, aac):
    d = Drive(S=S, alpha_ac=aac)
    return rocking_cube(alpha, d, np.array([d.q]), t, n_z=int(60 * aac))[:, 0]


print(f"(a) S scan, alpha = 3, alpha_ac = 60; zero at {Z:.4f}")
for dS in (-0.3, -0.1, -0.03, -0.01, 0.0, 0.01, 0.03, 0.1, 0.3):
    S = Z + dS
    I = centre(S, 3.0, 60.0)
    th = np.sin(np.pi * 3.0 * jv(1, S) * t) ** 2
    print(f"  S-Z={dS:+.2f} J1={jv(1,S):+.4f}  max I={I.max():.2e} (law {th.max():.2e})  "
          f"I(t=1)={I[-1]:.2e} (law {th[-1]:.2e})")

print("(b) residual at the zero: max over t of I(q, t)")
for alpha in (1.0, 2.0, 3.0):
    row = []
    for aac in (30.0, 60.0, 120.0):
        I = centre(Z, alpha, aac)
        row.append(f"aac={aac:>5}: {I.max():.2e}")
    print(f"  alpha={alpha}: " + "  ".join(row))

print("(c) minimum of max_t I over S near the zero (alpha = 3)")
for aac in (30.0, 60.0):
    Ss = Z + np.linspace(-0.06, 0.06, 25)
    m = [centre(S, 3.0, aac).max() for S in Ss]
    k = int(np.argmin(m))
    print(f"  aac={aac}: min at S-Z={Ss[k]-Z:+.4f}, value {m[k]:.2e}")
