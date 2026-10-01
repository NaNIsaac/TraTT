"""Entry 32: two co-propagating acoustic waves, the second launched tau later.

phi = S1 A1(s) sin(q1 s') + S2 A2(s) sin(q2 s'), s = t - zeta, both step envelopes,
A2 delayed by tau. Channel (n, m) at eta = n q1 + m q2 has coupling J_n(S1 A1) J_m(S2 A2)
(product of the two Jacobi-Anger expansions), so
    Theta_nm(t) = int_0^1 J_n(S1 A1(t - zeta)) J_m(S2 A2(t - zeta)) dzeta
and combination channels (n, m both != 0) can only grow where the two waves overlap.
Laue, alpha = 0.5; alpha_ac1 = 40, alpha_ac2 = 57 (incommensurate on the channel scale).
"""
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.core import rocking  # noqa: E402

ALPHA, A1, A2, S1, S2, TAU = 0.5, 40.0, 57.0, 1.2, 1.2, 0.4
q1, q2 = 2 * np.pi * A1, 2 * np.pi * A2


@dataclass(frozen=True)
class Two:
    t: float

    def __call__(self, zeta):
        s = self.t - np.asarray(zeta, dtype=np.float64)
        w1 = np.where(s >= 0, S1, 0.0)
        w2 = np.where(s >= TAU, S2, 0.0)
        return w1 * np.sin(-q1 * s) + w2 * np.sin(-q2 * s)


def theta(n, m, t, z=np.linspace(0, 1, 4001)):
    s = t - z
    return np.trapezoid(jv(n, S1 * (s >= 0)) * jv(m, S2 * (s >= TAU)), z)


chans = [(1, 0), (0, 1), (1, 1), (1, -1), (2, -1)]
eta = np.array([n * q1 + m * q2 for n, m in chans])
t = np.linspace(0.0, 1.0 + TAU + 0.2, 33) + 0.0037
nz = int(20 * (np.abs(eta).max() / (2 * np.pi) + (S1 + S2) * A2))
I = np.array([rocking(ALPHA, eta, phi=Two(x), n_z=nz) for x in t])
law = np.array([[np.sin(np.pi * ALPHA * abs(theta(n, m, x))) ** 2 for n, m in chans] for x in t])
print(f"n_z={nz}, law max err per channel:", dict(zip(chans, np.round(np.abs(I - law).max(0), 4))))
print("t     " + "  ".join(f"{str(c):>8s}" for c in chans))
for k in range(0, t.size, 3):
    print(f"{t[k]:.2f}  " + "  ".join(f"{v:8.4f}" for v in I[k]))
first = {c: float(t[np.argmax(I[:, j] > 0.2 * I[:, j].max())]) for j, c in enumerate(chans)}
print("time of first reaching 20 % of its maximum:", first)
