"""Entry 2: is the excess of p_1 over the collapse law the unmodulated rocking-curve tail?

Estimate: off Bragg |Dh|^2 ~ (2 pi alpha / eta)^2 sin^2(eta/2), fringe average 1/2,
integrated over the window (q/2, 3q/2) gives 8 (pi alpha)^2 / (3 q).
Test: (i) excess scales as alpha^2 / alpha_ac; (ii) subtracting the S = 0 window
integral removes it.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.core import rocking  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

S, n = 0.5, 1
t = np.array([0.25, 0.5, 0.75, 1.0])
for alpha in (1.0, 2.0):
    for aac in (10.0, 20.0, 40.0):
        d = Drive(S=S, alpha_ac=aac)
        q = d.q
        eta = np.linspace(0.5 * q, 1.5 * q, 2401)
        I = rocking_cube(alpha, d, eta, t, n_z=max(1200, int(30 * aac)))
        p = np.trapezoid(I, eta, axis=1)
        p0 = np.trapezoid(rocking(alpha, eta, n_z=1200), eta)
        an = alpha * abs(jv(n, S))
        th = np.array([np.pi**2 * an * quad(lambda x: jv(0, x), 0, 2 * np.pi * an * x)[0] for x in t])
        est = 8 * (np.pi * alpha) ** 2 / (3 * q)
        ex = p - th
        print(f"alpha={alpha} aac={aac:>4} est_tail={est:.3f} S0_window={p0:.3f} "
              f"excess={np.array2string(ex, precision=3)} "
              f"excess-S0={np.array2string(ex - p0, precision=3)} rel_after={np.max(np.abs(ex - p0) / th):.3f}")
