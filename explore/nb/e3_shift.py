"""Entry 3: where does the satellite peak sit, and why does it move?

Peak position from a parabola through the three highest samples of raw I on a
fine eta grid. Measured at t = 1 (crystal fully filled) and during transit.
Scalings distinguish candidates:
  tail slope pulling the peak   ~ alpha^2 / q^3 / (h t^2)   (artifact of reading)
  Floquet level shift           ~ (pi alpha J)^2 / q         (genuine resonance shift)
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.core import amplitude  # noqa: E402
from tratt.drive import Drive  # noqa: E402


def peakpos(eta, y):
    k = int(np.argmax(y))
    a, b, c = y[k - 1], y[k], y[k + 1]
    h = eta[1] - eta[0]
    return eta[k] + 0.5 * h * (a - c) / (a - 2 * b + c)


S = 1.5
for aac in (20.0, 40.0):
    d = Drive(S=S, alpha_ac=aac)
    q = d.q
    for alpha in (0.5, 1.0, 2.0):
        an = alpha * jv(1, S)
        out = []
        for t in (0.3, 0.6, 1.0):
            row = []
            for n in (1, -1):
                eta = np.linspace(n * q - 12, n * q + 12, 2401)
                A = amplitude(alpha, eta, phi=d.phi(t), n_z=max(1500, int(40 * aac)))
                A0 = amplitude(alpha, eta, n_z=max(1500, int(40 * aac)))
                row.append((peakpos(eta, np.abs(A) ** 2) - n * q, peakpos(eta, np.abs(A - A0) ** 2) - n * q))
            out.append(f"t={t}: raw(+1,-1)=({row[0][0]:+.3f},{row[1][0]:+.3f}) diff=({row[0][1]:+.3f},{row[1][1]:+.3f})")
        print(f"aac={aac} alpha={alpha} (pi a J0)^2/q={(np.pi*alpha*jv(0,S))**2/q:.4f} "
              f"(pi a)^2/q={(np.pi*alpha)**2/q:.4f} alpha_n={an:.2f}\n   " + "\n   ".join(out))
