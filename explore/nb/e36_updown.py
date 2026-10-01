"""Entry 36: same transducer amplitude, up vs down a triangular pulse — the whole I(eta).

Classification (user's three levels):
  1  only the relative satellite intensities differ   -> history-dependent intensity
  2  widths / splitting / fringes differ visibly      -> history-dependent line shape
  3  new peak positions or combination structure      -> genuinely new
Extra test: is each channel's line shape that of a *uniform* crystal with the same
accumulated coupling Theta_n (Laue slab, coupling pi alpha Theta_n, thickness 1)?
If yes, the depth history is compressed into one number per channel.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from memory_drive import HistoryDrive, theta, triangle  # noqa: E402
from tratt.core import required_nz, rocking  # noqa: E402

AAC, SMAX = 40.0, 3.0
NS = (-3, -2, -1, 0, 1, 2, 3)


def uniform_line(d, a):
    """Laue slab of thickness 1 and coupling a at detuning d."""
    s = np.sqrt(a * a + d * d / 4)
    return a * a * np.sin(s) ** 2 / s**2


def fwhm(x, y):
    k = int(np.argmax(y))
    h = y[k] / 2
    i, j = k, k
    while i > 0 and y[i] > h:
        i -= 1
    while j < y.size - 1 and y[j] > h:
        j += 1
    return x[j] - x[i]


def local_peaks(x, y, rel=0.05):
    return [round(float(x[i]), 2) for i in range(1, y.size - 1)
            if y[i] > y[i - 1] and y[i] >= y[i + 1] and y[i] > rel * y.max()]


# "partial": first run, where t_up < 1 (crystal not yet filled on the way up - confounded).
# Default: both moments fully filled, t_up = A T >= 1.
PARTIAL = [(0.5, 1.0, 0.6), (2.0, 1.0, 0.6), (2.0, 0.5, 0.5)]
FILLED = [(0.5, 2.0, 0.6), (2.0, 2.0, 0.6), (2.0, 2.0, 0.8), (2.0, 1.25, 0.8)]
CASES = PARTIAL if "partial" in sys.argv else FILLED
for alpha, T, A in CASES:
    d = HistoryDrive(S=SMAX, alpha_ac=AAC, env=triangle(T))
    q = d.q
    tu, td = A * T + 0.0037, (2 - A) * T + 0.0037
    eta = np.arange(-3.5 * q, 3.5 * q, 0.25)
    nz = required_nz(alpha, AAC, S=SMAX, eta_max=3.5 * q)
    Iu, Id = (rocking(alpha, eta, phi=d.phi(x), n_z=nz) for x in (tu, td))
    print(f"=== alpha={alpha} T={T}: transducer amplitude S A = {SMAX*A:.2f} at t_up={tu:.3f}, t_down={td:.3f}")
    print("  n | window p up/down | centre up/down | FWHM up/down (uniform-Theta FWHM) | "
          "shape diff | peaks within ±q/2 (up ; down)")
    for n in NS:
        m = np.abs(eta - n * q) < q / 2
        x = eta[m] - n * q
        pu, pd = np.trapezoid(Iu[m], x), np.trapezoid(Id[m], x)
        k0 = np.argmin(np.abs(x))
        wu, wd = fwhm(x, Iu[m]), fwhm(x, Id[m])
        thu, thd = abs(theta(d, n, [tu])[0]), abs(theta(d, n, [td])[0])
        wu0 = fwhm(x, uniform_line(x, np.pi * alpha * thu)) if thu > 1e-3 else np.nan
        wd0 = fwhm(x, uniform_line(x, np.pi * alpha * thd)) if thd > 1e-3 else np.nan
        diff = np.max(np.abs(Iu[m] / max(pu, 1e-12) - Id[m] / max(pd, 1e-12))) * max(wu, wd)
        print(f"  {n:+d} | {pu:7.3f}/{pd:7.3f} | {Iu[m][k0]:.3f}/{Id[m][k0]:.3f} | "
              f"{wu:6.2f}/{wd:6.2f} ({wu0:6.2f}/{wd0:6.2f}) | {diff:5.2f} | "
              f"{local_peaks(x, Iu[m])[:5]} ; {local_peaks(x, Id[m])[:5]}")
    gaps = [(float(eta[i]) / q) for i in range(1, eta.size - 1)
            if Iu[i] > 0.02 * Iu.max() and Iu[i] > Iu[i - 1] and Iu[i] >= Iu[i + 1]
            and min(abs(eta[i] / q - round(eta[i] / q)), 1) > 0.25]
    print(f"  peaks > 2 % of max outside ±q/4 of any channel (up): {np.round(gaps, 3)}")
