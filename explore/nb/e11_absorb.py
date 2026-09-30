"""Entry 11b: Borrmann absorption via complex alpha, centre intensity of satellite n.

Prediction (mapping unchanged, coupling now complex):
    I(n q, t) = |sin(pi alpha J_n(S) t)|^2 = [cosh(2 pi alpha'' |J_n| t) - cos(2 pi alpha' |J_n| t)] / 2
(common exp(-mu L) gauged out). Contrast of the oscillation ~ 1 / cosh(2 pi alpha'' |J_n| t).
Mapping to mu: alpha'' = alpha' * chi_h'' / chi_h', and for Si 220 chi_h''/chi_0'' ~ kappa.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.core import amplitude  # noqa: E402
from tratt.drive import Drive  # noqa: E402

t = np.linspace(0.0, 1.0, 101) + 0.0037
for S, a_re, a_im, n in [(1.5, 2.0, 0.05, 1), (1.5, 2.0, 0.15, 1), (3.0, 2.0, 0.1, 2), (1.5, 4.0, 0.1, 1)]:
    d = Drive(S=S, alpha_ac=20.0)
    alpha = complex(a_re, a_im)
    I = np.array([abs(amplitude(alpha, np.array([n * d.q]), phi=d.phi(x), n_z=1200)[0]) ** 2 for x in t])
    th = np.abs(np.sin(np.pi * alpha * abs(jv(n, S)) * t)) ** 2
    print(f"S={S} alpha={alpha} n={n}: max|I-th|={np.max(np.abs(I - th)):.3f}  "
          f"I(max) sim/th = {I.max():.3f}/{th.max():.3f}  min after 1st max sim/th = "
          f"{I[np.argmax(I):].min():.3f}/{th[np.argmax(th):].min():.3f}")
