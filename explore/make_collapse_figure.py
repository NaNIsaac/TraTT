"""Satellite formation collapses onto the static Laue thickness curve. Run: python explore/make_collapse_figure.py

Inside the filled depth zeta < t the grating is fully developed, and channel n is
a resonant two-beam problem with coupling pi*alpha*J_n(S); every other channel
is detuned by a multiple of q >> pi*alpha. The unfilled remainder is off-Bragg
at eta = n q. So the satellite is a Laue crystal of thickness t_tilde, and

    center     I(n q, t)  ~ sin^2(pi alpha_n t)
    integrated p_n(t)     ~ pi^2 alpha_n int_0^{2 pi alpha_n t} J_0(x) dx

with alpha_n = alpha |J_n(S)|. Both are plotted against alpha_n * t with no fit
parameter. Raw intensity is used throughout -- this is what a detector records.
The regime is resolved channels, 2 pi alpha << q (alpha_ac = 20 here).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

INK, INK2 = "#0b0b0b", "#52514e"
HUES = ["#2a78d6", "#eb6834", "#1baf7a", "#9b59b6", "#c0392b"]
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})

CASES = [(0.5, 2.0, 1), (1.5, 2.0, 1), (1.5, 4.0, 1), (2.0, 4.0, 2), (1.0, 6.0, 1)]  # (S, alpha, n)


def main(out="figures/satellite_collapse.png"):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    t = np.linspace(0.01, 1.0, 100)
    x = np.linspace(0, 2 * np.pi * 3.0, 4000)
    F = cumulative_trapezoid(jv(0, x), x, initial=0.0)
    for (S, alpha, n), c in zip(CASES, HUES):
        d = Drive(S=S, alpha_ac=20.0)
        q = d.q
        eta = np.linspace((n - 0.5) * q, (n + 0.5) * q, 1201)
        I = rocking_cube(alpha, d, eta, t, n_z=1000)
        an = alpha * abs(jv(n, S))
        k0 = np.argmin(np.abs(eta - n * q))
        p = np.trapezoid(I, eta, axis=1)
        lab = rf"$S$={S}, $\alpha$={alpha:g}, $n$={n}  ($\alpha_n$={an:.2f})"
        ax[0].plot(an * t, I[:, k0], "o", ms=2.5, color=c, label=lab)
        ax[1].plot(an * t, p / (np.pi**2 * an), "o", ms=2.5, color=c, label=lab)
    s = np.linspace(0, 3.0, 600)
    ax[0].plot(s, np.sin(np.pi * s) ** 2, color=INK, lw=1.2, label=r"$\sin^2(\pi\alpha_n\tilde t)$")
    ax[1].plot(x / (2 * np.pi), F, color=INK, lw=1.2,
               label=r"$\int_0^{2\pi\alpha_n\tilde t}J_0(x)\,dx$")
    ax[1].axhline(1.0, color=INK2, lw=0.8, ls=":")
    ax[0].set(xlabel=r"$\alpha_n\tilde t$", ylabel=r"$I(\eta=nq,\tilde t)$", xlim=(0, 2.3),
              title="a  Satellite centre: temporal Pendellösung")
    ax[1].set(xlabel=r"$\alpha_n\tilde t$", ylabel=r"$p_n/(\pi^2\alpha_n)$", xlim=(0, 2.3),
              title="b  Integrated satellite: collapse, no fit parameter")
    ax[0].legend(fontsize=7, frameon=False, loc="upper right")
    ax[1].legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle(r"Laue, step onset, $\alpha_{\rm ac}=20$, raw intensity, $\tilde t=t/t_{\rm tr}\leq 1$;"
                 r"  $\alpha_n=\alpha|J_n(S)|$", fontsize=10)
    Path(out).parent.mkdir(exist_ok=True)
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
