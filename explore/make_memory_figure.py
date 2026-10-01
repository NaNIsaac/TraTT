"""Path memory of satellite formation, four panels. Run: python explore/make_memory_figure.py

The crystal stores the last t_tr of the drive; satellite n reads
Theta_n(t) = int_0^1 J_n(S A(t - zeta)) sign_n dzeta (sin^2 in Laue, tanh^2 in Bragg).
Laue, alpha = 0.5 (no Pendellosung), alpha_ac = 40. Points: full TT; lines: the law.
a  transient-only satellite (S_f at the J_1 zero, ramp T = 2)
b  same transducer amplitude, different pattern: I_1 vs S A(t) for a triangular pulse
c  multiple revival: slow triangle to S = 10
d  phase-coded packet: pi flip in the middle cancels n = 1 while n = 2 peaks
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jn_zeros

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "explore" / "nb"))

from memory_drive import HistoryDrive, box, ramp, theta, triangle  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200, "lines.linewidth": 1.8,
})
ALPHA, AAC, OFF = 0.5, 40.0, 0.0037


def centre(d, ns, t):
    return rocking_cube(ALPHA, d, np.array([n * d.q for n in ns]), t, n_z=2400)


def law(d, n, t):
    return np.sin(np.pi * ALPHA * np.abs(theta(d, n, t))) ** 2


def pts(ax, t, y, c, m, lab=None):
    ax.plot(t, y, m, ms=3.5, color=c, mec="white", mew=0.5, ls="none", label=lab)


def main(out="figures/poster_memory.png"):
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.6), constrained_layout=True)

    # a: transient-only satellite
    d = HistoryDrive(S=jn_zeros(1, 1)[0], alpha_ac=AAC, env=ramp(2.0))
    t = np.linspace(0, 3.4, 69) + OFF
    tf = np.linspace(0, 3.4, 400)
    I = centre(d, (1, 2, 3), t)
    for j, (n, c, m) in enumerate(zip((1, 2, 3), (BLUE, ORANGE, AQUA), ("o", "s", "^"))):
        ax[0, 0].plot(tf, law(d, n, tf), color=c)
        pts(ax[0, 0], t, I[:, j], c, m, rf"$n={n}$")
    ax[0, 0].axvline(3.0, color=INK2, lw=0.8, ls="--")
    ax[0, 0].text(2.3, 0.6, "ramp has left\nthe crystal", fontsize=8, color=INK2)
    ax[0, 0].set(xlabel=r"$\tilde t$", ylabel=r"$I(nq,\tilde t)$", ylim=(-0.02, 0.72),
                 title=r"a  Transient-only satellite: $S_f=j_{1,1}$, $J_1(S_f)=0$")
    ax[0, 0].legend(fontsize=8, frameon=False, loc="upper right")

    # b: hysteresis loop, triangular pulse
    for T, c, m in ((1.0, BLUE, "o"), (4.0, ORANGE, "s")):
        d = HistoryDrive(S=3.0, alpha_ac=AAC, env=triangle(T))
        t = np.linspace(0, 2 * T, 81) + OFF
        I = centre(d, (1,), t)[:, 0]
        A = 3.0 * d.envelope(t)
        tf = np.linspace(0, 2 * T, 400)
        ax[0, 1].plot(3.0 * d.envelope(tf), law(d, 1, tf), color=c)
        pts(ax[0, 1], A, I, c, m, rf"rise = fall = {T:g} $t_{{\rm tr}}$")
        k = len(tf) * 3 // 4
        ax[0, 1].annotate("", xy=(3.0 * d.envelope(tf[k + 5]), law(d, 1, tf[k + 5:k + 6])[0]),
                          xytext=(3.0 * d.envelope(tf[k]), law(d, 1, tf[k:k + 1])[0]),
                          arrowprops=dict(arrowstyle="->", color=c, lw=1.2))
    ax[0, 1].set(xlabel=r"transducer amplitude $S\,A(\tilde t)$", ylabel=r"$I(q,\tilde t)$",
                 title="b  Same strain, different pattern: up vs down a triangular pulse")
    ax[0, 1].legend(fontsize=8, frameon=False, loc="upper left")

    # c: multiple revival
    d = HistoryDrive(S=10.0, alpha_ac=AAC, env=triangle(3.0))
    t = np.linspace(0, 7.0, 141) + OFF
    tf = np.linspace(0, 7.0, 700)
    ax[1, 0].plot(tf, law(d, 1, tf), color=BLUE)
    pts(ax[1, 0], t, centre(d, (1,), t)[:, 0], BLUE, "o", r"$n=1$")
    ax2 = ax[1, 0].twiny()  # same y quantity; top axis only labels the drive level
    ax2.set_xlim(ax[1, 0].get_xlim())
    ticks = [0, 1.15, 2.1, 3.0, 3.9, 4.85, 6.0]
    ax2.set_xticks(ticks, [f"{10*d.envelope(np.array(x)):.1f}" for x in ticks], fontsize=7)
    ax2.set_xlabel(r"transducer amplitude $S\,A$", fontsize=8, color=INK2)
    ax[1, 0].set(xlabel=r"$\tilde t$", ylabel=r"$I(q,\tilde t)$",
                 title=r"c  Multiple revival, slow pulse to $S$=10 (no Pendellösung, $\alpha$=0.5)")

    # d: phase-coded packet
    d = HistoryDrive(S=1.5, alpha_ac=AAC, env=box(0.8), flip=lambda s: np.where(s >= 0.4, 1.0, 0.0))
    t = np.linspace(0, 2.0, 81) + OFF
    tf = np.linspace(0, 2.0, 500)
    I = centre(d, (1, 2), t)
    for j, (n, c, m) in enumerate(zip((1, 2), (BLUE, ORANGE), ("o", "s"))):
        ax[1, 1].plot(tf, law(d, n, tf), color=c)
        pts(ax[1, 1], t, I[:, j], c, m, rf"$n={n}$")
    ax[1, 1].axvspan(0.8, 1.0, color=GRID, zorder=0)
    ax[1, 1].text(0.82, 0.125, "both halves\ninside", fontsize=8, color=INK2)
    ax[1, 1].set(xlabel=r"$\tilde t$", ylabel=r"$I(nq,\tilde t)$", ylim=(-0.005, 0.15),
                 title=r"d  Packet with a $\pi$ carrier flip: $n=1$ cancels, then revives")
    ax[1, 1].legend(fontsize=8, frameon=False, loc="upper right")

    for a in ax.flat:
        a.grid(color=GRID, lw=0.6)
        a.set_axisbelow(True)
    fig.suptitle(r"Formation remembers the path: satellite $n$ reads $\Theta_n(\tilde t)=\int_0^1 "
                 r"J_n(S A(\tilde t-\zeta))\,d\zeta$ (Laue, $\alpha$=0.5, $\alpha_{\rm ac}$=40;"
                 r" points full TT, lines closed form)", fontsize=10)
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
