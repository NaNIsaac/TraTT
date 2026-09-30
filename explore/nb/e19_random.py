"""Entry 19: random-condition test of the centre law (guards against special points).

30 conditions drawn with a fixed seed: S in [0.3, 3.5], alpha in [0.3, 4], alpha_ac
non-integer in [20, 45], n in {1, 2, -1}, rise in {0, U[0, 0.6]}. Law: accumulated
coupling (Entry 5), which is the step law when rise = 0. Report the error against
alpha / alpha_ac.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

rng = np.random.default_rng(20260930)
z = np.linspace(0.0, 1.0, 2001)
rows = []
for i in range(30):
    S = rng.uniform(0.3, 3.5)
    alpha = rng.uniform(0.3, 4.0)
    aac = rng.uniform(20.0, 45.0)
    n = int(rng.choice([1, 2, -1]))
    rise = 0.0 if rng.random() < 0.5 else rng.uniform(0.0, 0.6)
    d = Drive(S=S, alpha_ac=aac, rise=rise)
    t = np.sort(rng.uniform(0.0, 1.0 + rise, 40))
    I = rocking_cube(alpha, d, np.array([n * d.q]), t, n_z=int(60 * aac))[:, 0]
    Th = np.array([abs(np.trapezoid(jv(n, S * d.envelope(x - z)), z)) for x in t])
    err = np.max(np.abs(I - np.sin(np.pi * alpha * Th) ** 2))
    rows.append((alpha / aac, err))
    print(f"{i:2d} S={S:.2f} alpha={alpha:.2f} aac={aac:.2f} n={n:+d} rise={rise:.2f} "
          f"ratio={alpha/aac:.3f} max err={err:.3f}")
r = np.array(rows)
for lo, hi in [(0, 0.05), (0.05, 0.1), (0.1, 0.2)]:
    m = (r[:, 0] > lo) & (r[:, 0] <= hi)
    if m.any():
        print(f"ratio in ({lo}, {hi}]: {m.sum()} cases, worst {r[m, 1].max():.3f}, median {np.median(r[m, 1]):.3f}")
