"""Entry 13: recover the acoustic amplitude S from temporal periods alone.

Simulated "measurement": centre intensity of n = 1 and n = 2 vs t in [0, 1], plus
multiplicative noise and an unknown overall scale (flux/gain). Fit each trace
with c * sin^2(pi f t) (f free, c free), then solve |J_2(S)/J_1(S)| = f_2/f_1 for S
on the first branch (S < 3.83, first zero of J_1). alpha is NOT used.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, curve_fit
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

rng = np.random.default_rng(0)
t = np.linspace(0.0, 1.0, 81) + 0.0037


def model(t, c, f):
    return c * np.sin(np.pi * f * t) ** 2


for S_true, alpha in [(1.0, 3.0), (1.5, 3.0), (2.0, 3.0), (2.5, 4.0)]:
    d = Drive(S=S_true, alpha_ac=20.0)
    I = rocking_cube(alpha, d, np.array([d.q, 2 * d.q]), t, n_z=1500)
    fs = []
    for j, n in enumerate((1, 2)):
        y = 0.37 * I[:, j] * (1 + 0.03 * rng.standard_normal(t.size))  # unknown gain, 3 % noise
        f0 = alpha * abs(jv(n, S_true)) * rng.uniform(0.8, 1.2)       # rough starting guess only
        best = min((curve_fit(model, t, y, p0=(y.max(), f), maxfev=20000) for f in
                    np.linspace(0.3 * f0, 1.7 * f0, 15)),
                   key=lambda r: np.sum((model(t, *r[0]) - y) ** 2))
        fs.append(abs(best[0][1]))
    r = fs[1] / fs[0]
    S_fit = brentq(lambda s: jv(2, s) / jv(1, s) - r, 0.05, 3.8)
    print(f"S_true={S_true} alpha={alpha}: f1={fs[0]:.3f} (th {alpha*abs(jv(1,S_true)):.3f}) "
          f"f2={fs[1]:.3f} (th {alpha*abs(jv(2,S_true)):.3f})  S_recovered={S_fit:.3f}  "
          f"error={100*(S_fit/S_true-1):+.1f} %")
