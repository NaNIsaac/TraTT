"""Entry 1: formation order of channels, and the two-layer law of the n = 0 channel.

Predictions (no fit parameter):
  n != 0 centre   I(nq, t) = sin^2(pi alpha |J_n(S)| t)   first max at t* = 1 / (2 alpha_n)
  n  = 0 centre   I(0, t)  = sin^2(pi alpha [1 - (1 - J_0(S)) t])
The n = 0 law: filled layer (coupling pi alpha J_0) and unfilled layer (pi alpha)
are both on Bragg at eta = 0, so their Laue rotations add.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

t = np.linspace(0.0, 1.0, 201)
for S, alpha in [(0.5, 2.0), (1.0, 2.0), (3.0, 2.0), (3.0, 1.0)]:
    d = Drive(S=S, alpha_ac=20.0)
    q = d.q
    ns = [0, 1, 2, 3]
    eta = np.array([n * q for n in ns])
    I = rocking_cube(alpha, d, eta, t, n_z=1200)
    print(f"S={S} alpha={alpha}")
    for j, n in enumerate(ns):
        J = jv(n, S)
        th = (np.sin(np.pi * alpha * (1 - (1 - J) * t)) ** 2 if n == 0
              else np.sin(np.pi * alpha * abs(J) * t) ** 2)
        err = np.max(np.abs(I[:, j] - th))
        # first local maximum of the simulated centre intensity (n != 0)
        y = I[:, j]
        k = int(np.argmax(y)) if n and 1 / (2 * alpha * abs(J)) < 1 else None
        tstar = t[k] if k is not None else np.nan
        pred = 1 / (2 * alpha * abs(J)) if n else np.nan
        print(f"  n={n} J_n={J:+.3f} max|I-th|={err:.3f}  first max sim t*={tstar:.3f} pred={pred:.3f}")
