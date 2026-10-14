"""Entry 65: depth-step convergence of the Si 400 TT check (single mode at eta = W).

RK4 at 24 steps per cycle of eta carries a phase error per step of order (h eta)^5 / 120; over
tens of thousands of cycles it accumulates. Run the same point at 24, 48 and 96 steps per cycle.
Usage: python e65_conv.py g x [x ...]
"""
import sys

import numpy as np

from e65_si400 import H, root
from tratt.core import solve_bragg

g = float(sys.argv[1])
alpha, kappa, nt = 6.0 / (2 * np.pi * g), 6.0, 32
S = float(__import__("os").environ.get("E65_S", "0.02"))   # drive amplitude (linear regime test)
for x in map(float, sys.argv[2:]):
    W = x * 2 * np.pi * alpha
    t = 2 * np.pi * np.arange(nt) / (W * nt)
    F_cf = H(x, g, alpha=alpha, S=S)[0] * abs(root(alpha, W, kappa)) * np.pi * alpha * S / (2 * kappa)
    for per in map(int, __import__("os").environ.get("E65_PER", "24,48,96").split(",")):
        nz = int(per * (W / (2 * np.pi) + alpha)) + 256
        R = np.abs(solve_bragg(alpha, np.array([W]), phi=lambda z: S * np.cos(W * (t[:, None] - z)),
                               n_z=nz, kappa=kappa)) ** 2
        F = abs(np.mean((R - R.mean(axis=0)) * np.exp(1j * W * t[:, None]), axis=0)[0])
        print(f"g={g} x={x:g} per_cycle={per:3d} n_z={nz:7d}: TT / closed form = {F / F_cf:.4f}", flush=True)
