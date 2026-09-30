"""Entry 3b: scaling of the raw satellite peak shift at t = 0.3 (off the period grid)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.core import amplitude  # noqa: E402
from tratt.drive import Drive  # noqa: E402


def peakpos(eta, y):
    k = int(np.argmax(y))
    a, b, c = y[k - 1], y[k], y[k + 1]
    return eta[k] + 0.5 * (eta[1] - eta[0]) * (a - c) / (a - 2 * b + c)


t = 0.3137  # not a whole number of acoustic periods for any alpha_ac used
for S in (1.5, -1.5):
    for aac in (20.0, 40.0, 80.0):
        d = Drive(S=S, alpha_ac=aac)
        q = d.q
        row = []
        for alpha in (0.01, 0.1, 0.5, 1.0):
            nz = max(1500, int(40 * aac))
            s = []
            for n in (1, -1):
                eta = np.linspace(n * q - 15, n * q + 15, 3001)
                A = amplitude(alpha, eta, phi=d.phi(t), n_z=nz)
                A0 = amplitude(alpha, eta, n_z=nz)
                s.append((peakpos(eta, abs(A) ** 2) - n * q, peakpos(eta, abs(A - A0) ** 2) - n * q))
            row.append(f"a={alpha}: raw {s[0][0]:+.3f}/{s[1][0]:+.3f} diff {s[0][1]:+.3f}/{s[1][1]:+.3f}")
        print(f"S={S:+} aac={aac}: " + " | ".join(row))
