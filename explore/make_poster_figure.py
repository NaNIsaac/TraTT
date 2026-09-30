"""Poster summary figure: the filled-depth mapping in four panels. Run: python explore/make_poster_figure.py

a  Laue centre intensity collapses on sin^2(pi alpha_n t)        (C-1.10, notebook Entry 1)
b  Bragg centre intensity collapses on tanh^2(pi alpha_n t)      (Entry 9)
c  Formation order follows the Bessel ladder, reorderable by S   (Entry 1)
d  Finite rise: sin^2 of the accumulated coupling                 (Entry 5)

All curves are closed forms with no fit parameter; points are full Takagi-Taupin.
Resolved channels (alpha_ac = 20), step onset unless stated, times sampled off
the acoustic period grid. Series carry marker shape as well as colour.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
HUES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # reference categorical order, validated
MARK = ["o", "s", "^", "D"]
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200, "lines.linewidth": 2,
})
AAC = 20.0
OFF = 0.0037  # keeps samples off whole acoustic periods (notebook Entry 2)
CASES = [(0.5, 2.0, 1), (1.5, 2.0, 1), (2.0, 4.0, 2), (1.5, 4.0, 1)]  # (S, alpha, n)


def centre(alpha, S, n, t, geometry="laue", rise=0.0):
    d = Drive(S=S, alpha_ac=AAC, rise=rise)
    return rocking_cube(alpha, d, np.array([n * d.q]), t, geometry=geometry, n_z=1200)[:, 0]


def collapse(ax, geometry, law, lawlabel, title):
    t = np.linspace(0.0, 1.0, 61) + OFF
    for (S, alpha, n), c, m in zip(CASES, HUES, MARK):
        an = alpha * abs(jv(n, S))
        ax.plot(an * t, centre(alpha, S, n, t, geometry), m, ms=4, mfc=c, mec="white", mew=0.6,
                color=c, ls="none", label=rf"$S$={S:g}, $\alpha$={alpha:g}, $n$={n}")
    s = np.linspace(0, 2.3, 500)
    ax.plot(s, law(s), color=INK, lw=1.4, label=lawlabel, zorder=0)
    ax.set(xlim=(0, 2.3), ylim=(-0.03, 1.08), xlabel=r"$\alpha_n\tilde t$,  $\alpha_n=\alpha|J_n(S)|$",
           ylabel=r"$I(\eta=nq,\,\tilde t)$", title=title)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right" if geometry == "bragg" else "upper right")


def order_panel(ax):
    alpha = 2.0
    Ss = np.linspace(0.3, 4.5, 400)
    t = np.linspace(0.0, 1.0, 201) + OFF
    for n, c, m in zip((1, 2, 3), HUES, MARK):
        tstar = 1 / (2 * alpha * np.abs(jv(n, Ss)))
        ax.plot(Ss, tstar, color=c, label=rf"$n={n}$: $1/(2\alpha|J_{n}(S)|)$")
        pts = []
        for S in (1.0, 2.0, 3.0, 4.0):
            if 1 / (2 * alpha * abs(jv(n, S))) < 0.95:
                pts.append((S, t[np.argmax(centre(alpha, S, n, t))]))
        if pts:
            ax.plot(*zip(*pts), m, ms=7, mfc=c, mec="white", mew=0.8, ls="none")
    ax.axhspan(1.0, 2.0, color=GRID, zorder=0)
    ax.text(0.35, 1.04, "after transit: no transient maximum", fontsize=8, color=INK2)
    ax.axvline(3.0, color=INK2, lw=0.8, ls=":")
    ax.annotate("S = 3: order n = 2, 1, 3", (3.0, 0.5), xytext=(3.15, 0.3), fontsize=8, color=INK2,
                arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))
    ax.set(xlim=(0.3, 4.5), ylim=(0.2, 1.2), xlabel=r"modulation index $S$",
           ylabel=r"first maximum $\tilde t^{*}$",
           title=r"c  Formation order follows $|J_n(S)|$ ($\alpha$=2; points: simulation)")
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")


def rise_panel(ax):
    S, alpha, rise = 1.5, 2.0, 0.5
    d = Drive(S=S, alpha_ac=AAC, rise=rise)
    t = np.linspace(0.0, 1.0 + rise, 91) + OFF
    z = np.linspace(0.0, 1.0, 2001)
    I = rocking_cube(alpha, d, np.array([0.0, d.q, 2 * d.q]), t, n_z=1200)
    tf = np.linspace(0.0, 1.0 + rise, 400)
    for j, (n, c, m) in enumerate(zip((0, 1, 2), HUES, MARK)):
        Th = np.array([abs(np.trapezoid(jv(n, S * d.envelope(x - z)), z)) for x in tf])
        ax.plot(tf, np.sin(np.pi * alpha * Th) ** 2, color=c, lw=1.6)
        ax.plot(t, I[:, j], m, ms=3.5, mfc=c, mec="white", mew=0.5, ls="none", label=rf"$n={n}$")
    Th1 = np.array([abs(np.trapezoid(jv(1, S * (x - z >= 0)), z)) for x in tf])
    ax.plot(tf, np.sin(np.pi * alpha * Th1) ** 2, color=INK2, lw=1.2, ls="--", label="$n=1$, step-onset law")
    ax.axvline(1.0 + rise, color=INK2, lw=0.8, ls=":")
    ax.set(xlim=(0, 1.0 + rise), ylim=(-0.03, 1.08), xlabel=r"$\tilde t=t/t_{\rm tr}$",
           ylabel=r"$I(\eta=nq,\,\tilde t)$",
           title=r"d  Rise 0.5: $\sin^2(\pi\alpha\,|\int_0^1 J_n(S\,A(\tilde t-\zeta))\,d\zeta|)$")
    ax.legend(fontsize=7.5, frameon=False, loc="upper left", ncol=2)


def main(out="figures/poster_mapping.png"):
    fig, ax = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    collapse(ax[0, 0], "laue", lambda s: np.sin(np.pi * s) ** 2, r"$\sin^2(\pi\alpha_n\tilde t)$",
             "a  Laue: temporal Pendellösung of each satellite")
    collapse(ax[0, 1], "bragg", lambda s: np.tanh(np.pi * s) ** 2, r"$\tanh^2(\pi\alpha_n\tilde t)$",
             "b  Bragg: saturation, formation time $\\approx 0.95\\,t_{\\rm tr}/\\alpha_n$")
    order_panel(ax[1, 0])
    rise_panel(ax[1, 1])
    for a in ax.flat:
        a.grid(color=GRID, lw=0.6)
        a.set_axisbelow(True)
    fig.suptitle("Satellite n forms as a crystal of the filled depth with coupling "
                 r"$\pi\alpha|J_n(S)|$ — curves are closed forms, no fit;"
                 r" points are full Takagi–Taupin ($\alpha_{\rm ac}$=20)", fontsize=10)
    fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
