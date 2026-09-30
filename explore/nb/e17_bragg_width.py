"""Entry 17: Bragg satellite width during formation.

Non-absorbing Bragg slab, coupling a = pi alpha_n, thickness t, u = (eta - n q)/2:
  |u| < a:  R = a^2 sinh^2(s t) / (s^2 + a^2 sinh^2(s t)),  s = sqrt(a^2 - u^2)
  |u| > a:  R = a^2 sin^2(s t)  / (s^2 + a^2 sin^2(s t)),   s = sqrt(u^2 - a^2)
t * FWHM -> 5.566 for a t -> 0; FWHM -> Darwin-like width ~ const * alpha_n for a t >> 1.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402
from tratt.observe import peak_fwhm  # noqa: E402


def bragg_slab(d, a, t):
    u = d / 2
    out = np.empty_like(u)
    inb = np.abs(u) < a
    s = np.sqrt(a * a - u[inb] ** 2)
    sh = np.sinh(s * t) ** 2
    out[inb] = a * a * sh / (s**2 + a * a * sh)
    s = np.sqrt(u[~inb] ** 2 - a * a)
    sn = np.sin(s * t) ** 2
    out[~inb] = a * a * sn / (s**2 + a * a * sn)
    return out


t = np.linspace(0.05, 1.0, 39) + 0.0037
for S, alpha, n in [(1.5, 1.0, 1), (1.5, 2.0, 1), (1.5, 4.0, 1)]:
    d = Drive(S=S, alpha_ac=20.0)
    q = d.q
    eta = np.linspace((n - 0.5) * q, (n + 0.5) * q, 4001)
    I = rocking_cube(alpha, d, eta, t, geometry="bragg", n_z=1500)
    an = alpha * abs(jv(n, S))
    dd = eta - n * q
    w_sim = np.array([peak_fwhm(eta, I[k], n * q, q / 2)[0] for k in range(t.size)])
    w_th = np.array([peak_fwhm(dd, bragg_slab(dd, np.pi * an, x), 0.0, q / 2)[0] for x in t])
    m = an * t > 0.3
    print(f"S={S} alpha={alpha} alpha_n={an:.2f}: max rel dev (alpha_n t>0.3) = "
          f"{np.nanmax(np.abs(w_sim[m] / w_th[m] - 1)):.3f};  FWHM/alpha_n at t=1: "
          f"sim {w_sim[-1]/an:.2f} th {w_th[-1]/an:.2f}")
    for k in range(0, t.size, 6):
        print(f"   t={t[k]:.2f} alpha_n t={an*t[k]:.2f}  FWHM sim {w_sim[k]:6.2f}  th {w_th[k]:6.2f}  "
              f"5.566/t {5.566/t[k]:6.2f}")
