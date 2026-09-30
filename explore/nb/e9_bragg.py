"""Entry 9: does the filled-depth mapping carry over to Bragg geometry?

Bragg, entrance face zeta = 0, acoustic front also enters at zeta = 0. At eta = n q
the filled layer [0, t] is a Bragg slab of coupling pi alpha_n; the unfilled layer
behind it is off-Bragg. Centre reflectivity of a non-absorbing Bragg slab of
thickness t at exact Bragg: tanh^2(pi alpha_n t).
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

t = np.linspace(0.0, 1.0, 101) + 0.0037
for S, alpha, n in [(1.5, 0.5, 1), (1.5, 1.0, 1), (1.5, 2.0, 1), (3.0, 2.0, 2)]:
    d = Drive(S=S, alpha_ac=20.0)
    q = d.q
    I = rocking_cube(alpha, d, np.array([n * q]), t, geometry="bragg", n_z=1500)[:, 0]
    an = alpha * abs(jv(n, S))
    th = np.tanh(np.pi * an * np.minimum(t, 1)) ** 2
    k = [10, 30, 60, 99]
    print(f"S={S} alpha={alpha} n={n} alpha_n={an:.2f} max|I-tanh2|={np.max(np.abs(I - th)):.3f}  "
          + " ".join(f"t={t[i]:.2f}:{I[i]:.3f}/{th[i]:.3f}" for i in k))
