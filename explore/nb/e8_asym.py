"""Entry 8: is the transient +-n asymmetry (C-1.7c) robust or a finite-alpha_ac interference?

A_1 = (p_{+1} - p_{-1}) / (p_{+1} + p_{-1}), window integrals of raw I.
The collapse law is symmetric in n (|J_{-n}| = |J_n|), so any asymmetry is beyond it.
Scan alpha_ac; also report the intensity-weighted rms over the transit.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402
from tratt.observe import channel_stats  # noqa: E402

S = 0.5
for alpha in (0.05, 1.3):
    for aac in (5.0, 10.0, 20.0, 40.0):
        d = Drive(S=S, alpha_ac=aac)
        q = d.q
        eta = np.linspace(-1.5 * q, 1.5 * q, 3001)
        t = np.linspace(0.005, 0.995, 199)
        I = rocking_cube(alpha, d, eta, t, n_z=max(800, int(40 * aac)))
        pp, pm = channel_stats(eta, I, d, 1).p, channel_stats(eta, I, d, -1).p
        A = (pp - pm) / (pp + pm)
        w = pp + pm
        rms_w = np.sqrt(np.sum(w * A**2) / np.sum(w))
        print(f"alpha={alpha} aac={aac:>4}: max|A|={np.max(np.abs(A)):.3f} at t={t[np.argmax(np.abs(A))]:.3f}  "
              f"weighted rms={rms_w:.4f}  rms(t>0.3)={np.sqrt(np.mean(A[t > 0.3]**2)):.4f}")
