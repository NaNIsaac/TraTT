"""Entry 41: validation of the strongest effects of Entries 38-40.

For each candidate: integrated satellite and peak with (base grid) / (2 n_z) / (eta step / 2)
/ (window x 1.5) / (wall smoothed over one and three acoustic wavelengths, where applicable),
plus the kinematic integral and the uniform-q control.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from qmod import PhaseDrive, kinematic, phase_wall, q_wall, spatial_chirp, uniform  # noqa: E402
from tratt.core import required_nz, rocking  # noqa: E402

AAC, S = 40.0, 1.5
Q0 = 2 * np.pi * AAC
T = 1.5037
W1 = 1.0 / AAC


def measure(alpha, theta_fn, qmax, half=0.45, step=0.1, nzmul=1):
    eta = np.arange(Q0 * (1 - half), Q0 * (1 + half) + (qmax - Q0), step)
    nz = nzmul * required_nz(alpha, qmax / (2 * np.pi), S=S, eta_max=eta.max())
    d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=theta_fn)
    I = rocking(alpha, eta, phi=d.phi(T), n_z=nz) - rocking(alpha, eta, n_z=nz)
    return np.trapezoid(I, eta), I.max(), d, eta


CANDIDATES = [
    ("chirp beta=0.4, alpha=2", 2.0, lambda w: spatial_chirp(Q0, 0.4), Q0 * 1.2, False),
    ("chirp beta=0.1, alpha=4", 4.0, lambda w: spatial_chirp(Q0, 0.1), Q0 * 1.05, False),
    ("pi wall zb=0.3, alpha=2", 2.0, lambda w: phase_wall(Q0, np.pi, T - 0.3, w), Q0, True),
    ("pi wall zb=0.5, alpha=2", 2.0, lambda w: phase_wall(Q0, np.pi, T - 0.5, w), Q0, True),
    ("q wall dc=2 zb=0.5, alpha=2", 2.0, lambda w: q_wall(Q0, Q0 + 4 * np.pi, T - 0.5, w), Q0 + 4 * np.pi, True),
]

if __name__ == "__main__":
    for alpha in (2.0, 4.0):
        p_u, pk_u, _, _ = measure(alpha, uniform(Q0), Q0)
        print(f"uniform control alpha={alpha}: integral {p_u:.3f}, peak {pk_u:.3f}")
    for name, alpha, mk, qmax, wall in CANDIDATES:
        base = measure(alpha, mk(0.0), qmax)
        res = [("base", base), ("2 n_z", measure(alpha, mk(0.0), qmax, nzmul=2)),
               ("eta/2", measure(alpha, mk(0.0), qmax, step=0.05)),
               ("window x1.5", measure(alpha, mk(0.0), qmax, half=0.6))]
        if wall:
            res += [("smooth 1 wl", measure(alpha, mk(W1), qmax)), ("smooth 3 wl", measure(alpha, mk(3 * W1), qmax))]
        d, eta = base[2], base[3]
        kin = np.trapezoid(kinematic(alpha, d, eta, T) - kinematic(alpha, PhaseDrive(S=0.0, alpha_ac=AAC), eta, T), eta)
        print(f"{name}: " + " | ".join(f"{k}: int {v[0]:.3f} peak {v[1]:.3f}" for k, v in res)
              + f" | kinematic int {kin:.3f}")
