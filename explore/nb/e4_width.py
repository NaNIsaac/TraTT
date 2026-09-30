"""Entry 4: spectral width of a forming satellite.

Theory: satellite n is a Laue slab of thickness t and coupling a = pi alpha_n,
    f(d) = a^2 sin^2(t sqrt(a^2 + d^2/4)) / (a^2 + d^2/4),   d = eta - n q
so t * FWHM = G(alpha_n t) is universal. G(0) = 5.566 (finite-depth sinc^2).
Measured on raw I with the same FWHM routine used for the theory curve.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402
from tratt.observe import peak_fwhm  # noqa: E402


def laue(d, a, t):
    s = np.sqrt(a * a + d * d / 4)
    return a * a * np.sin(t * s) ** 2 / s**2


def G(x, halfwin=200.0):
    """t * FWHM of the Laue curve at alpha_n t = x (t = 1 without loss)."""
    d = np.linspace(-halfwin, halfwin, 40001)
    w, _, _ = peak_fwhm(d, laue(d, np.pi * x, 1.0), 0.0, halfwin)
    return w


if __name__ == "__main__":
    xs = np.linspace(0.02, 1.4, 70)
    Gs = np.array([G(x) for x in xs])
    print("theory t*W:", " ".join(f"{x:.2f}:{g:.2f}" for x, g in zip(xs[::7], Gs[::7])))
    t = np.linspace(0.1, 1.0, 46) + 0.0037  # off the acoustic period grid
    for S, alpha, n in [(1.5, 1.0, 1), (1.5, 2.0, 1), (3.0, 2.0, 2), (0.5, 4.0, 1)]:
        d = Drive(S=S, alpha_ac=20.0)
        q = d.q
        eta = np.linspace((n - 0.5) * q, (n + 0.5) * q, 4001)
        I = rocking_cube(alpha, d, eta, t, n_z=1200)
        an = alpha * abs(jv(n, S))
        tw = np.array([peak_fwhm(eta, I[k], n * q, q / 2)[0] * t[k] for k in range(t.size)])
        th = np.interp(an * t, xs, Gs)
        m = an * t < 0.9
        print(f"S={S} a={alpha} n={n} alpha_n={an:.2f}  max rel dev (alpha_n t<0.9) = "
              f"{np.nanmax(np.abs(tw[m] / th[m] - 1)):.3f}")
        for k in range(0, t.size, 6):
            print(f"   alpha_n t={an*t[k]:.2f}  sim t*W={tw[k]:.2f}  th={th[k]:.2f}")
