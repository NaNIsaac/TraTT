"""Up vs down a triangular pulse at the same transducer amplitude: the whole rocking curve.
Run: python explore/make_updown_figure.py

Laue, alpha = 2, alpha_ac = 40, S_max = 3, rise = fall = 2 t_tr; S A = 2.4 at t = 1.6 (up)
and t = 2.4 (down); the crystal is fully filled at both moments (notebook Entry 36).
a  I(eta) over seven channels, up and down overlaid.
b  area-normalised line shapes of n = 2 and n = 3, with the uniform-crystal line shape
   that has the same accumulated coupling Theta_n (dotted).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "explore" / "nb"))

from memory_drive import HistoryDrive, theta, triangle  # noqa: E402
from tratt.core import required_nz, rocking  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE = "#2a78d6", "#eb6834"
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})


def uniform_line(d, a):
    s = np.sqrt(a * a + d * d / 4)
    return a * a * np.sin(s) ** 2 / s**2


def main(out="figures/updown_rocking.png"):
    alpha, T, A, SMAX = 2.0, 2.0, 0.8, 3.0
    d = HistoryDrive(S=SMAX, alpha_ac=40.0, env=triangle(T))
    q = d.q
    tu, td = A * T + 0.0037, (2 - A) * T + 0.0037
    eta = np.arange(-3.5 * q, 3.5 * q, 0.25)
    nz = required_nz(alpha, 40.0, S=SMAX, eta_max=3.5 * q)
    Iu, Id = (rocking(alpha, eta, phi=d.phi(x), n_z=nz) for x in (tu, td))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True,
                           gridspec_kw=dict(width_ratios=[2.2, 1, 1]))
    ax[0].plot(eta / q, Iu, color=BLUE, lw=1.0, label=rf"way up ($\tilde t$={tu:.2f})")
    ax[0].plot(eta / q, Id, color=ORANGE, lw=1.0, label=rf"way down ($\tilde t$={td:.2f})")
    ax[0].set(xlabel=r"$\eta/q$ (channel index)", ylabel=r"$I(\eta)$", xlim=(-3.5, 3.5),
              title=rf"a  Same transducer amplitude $SA$={SMAX*A:.1f}, two histories")
    ax[0].legend(fontsize=8, frameon=False, loc="upper right")
    for i, n in enumerate((2, 3)):
        m = np.abs(eta - n * q) < 30
        x = eta[m] - n * q
        for I, c, lab, t in ((Iu, BLUE, "up", tu), (Id, ORANGE, "down", td)):
            p = np.trapezoid(I[m], x)
            ax[1 + i].plot(x, I[m] / p, color=c, lw=1.6, label=f"{lab}  (area {p:.2f})")
            th = abs(theta(d, n, [t])[0])
            u = uniform_line(x, np.pi * alpha * th)
            ax[1 + i].plot(x, u / np.trapezoid(u, x), color=c, lw=1.0, ls=":")
        ax[1 + i].set(xlabel=rf"$\eta-{n}q$", ylabel="area-normalised I",
                      title=rf"{'bc'[i]}  $n={n}$; dotted: uniform, same $\Theta_{n}$")
        ax[1 + i].legend(fontsize=8, frameon=False, loc="upper right")
    for a in ax:
        a.grid(color=GRID, lw=0.6)
        a.set_axisbelow(True)
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
