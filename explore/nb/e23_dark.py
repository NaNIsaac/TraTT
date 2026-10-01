"""Entry 23: where does the intensity go when the satellite centre goes dark?

Track, vs t: centre I(nq), satellite-window integral, all-eta total, and the two side
lobes eta_+-(t) - n q. Theory (satellite = Laue slab, coupling a = pi alpha_n, thickness t,
detuning d = eta - n q):  f(d) = a^2 sin^2(t sqrt(a^2 + d^2/4)) / (a^2 + d^2/4)
so the lobe positions obey (eta_+- - n q) * t = G_+-(alpha_n t)  (universal).
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402


def laue_line(d, a, t):
    s = np.sqrt(a * a + d * d / 4)
    return a * a * np.sin(t * s) ** 2 / s**2


def lobes(dd, y, half):
    """Positions of the highest local maximum on each side of d = 0 (|d| < half)."""
    out = []
    for side in (-1, 1):
        m = (side * dd > 0.2) & (np.abs(dd) < half)
        e, v = dd[m], y[m]
        pk = [i for i in range(1, v.size - 1) if v[i] >= v[i - 1] and v[i] > v[i + 1]]
        out.append(e[max(pk, key=lambda i: v[i])] if pk else np.nan)
    return out


def run(S, alpha, aac, n, label):
    d = Drive(S=S, alpha_ac=aac)
    q = d.q
    an = alpha * abs(jv(n, S))
    eta = np.arange(-2.5 * q, 2.5 * q, 0.1)
    t = np.linspace(0.05, 1.0, 39) + 0.0037
    I = rocking_cube(alpha, d, eta, t, n_z=int(60 * aac))
    dd = eta - n * q
    win = np.abs(dd) < q / 2
    k0 = np.argmin(np.abs(dd))
    total = np.trapezoid(I, eta, axis=1)
    p = np.trapezoid(I[:, win], eta[win], axis=1)
    print(f"--- {label}: S={S} alpha={alpha} aac={aac} n={n} alpha_n={an:.3f} "
          f"(dark at t={1/an:.3f})")
    print(" alpha_n t | centre  th    | window p/(pi^2 a_n)  th   | total  | lobes (eta-nq)*t sim  th")
    for k in range(0, t.size, 2):
        x = an * t[k]
        thp = quad(lambda u: jv(0, u), 0, 2 * np.pi * x)[0]
        lo = lobes(dd, I[k], q / 2)
        dth = np.linspace(-q / 2, q / 2, 20001)
        lth = lobes(dth, laue_line(dth, np.pi * an, t[k]), q / 2)
        print(f"  {x:6.3f}  | {I[k, k0]:.3f} {np.sin(np.pi*x)**2:.3f} | "
              f"{p[k]/(np.pi**2*an):6.3f}            {thp:6.3f} | {total[k]:7.3f} | "
              f"{lo[0]*t[k]:+6.2f} {lo[1]*t[k]:+6.2f}   {lth[0]*t[k]:+6.2f} {lth[1]*t[k]:+6.2f}")


if __name__ == "__main__":
    run(1.5, 3.0, 60.0, 1, "map case")
    run(1.5, 2.0, 40.0, 1, "second case")
