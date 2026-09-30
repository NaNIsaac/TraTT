"""Poster figures. Run: python explore/make_figures.py

One condition in, four panels out. Colors follow the categorical slots 1-3
(blue / orange / aqua), which validate on all pairs for both normal and CVD
vision; +n and -n share a hue and are separated by line style, so channel
identity is never carried by color alone. Every line is directly labeled, which
is also the relief the aqua slot needs for its sub-3:1 contrast.

Intensity maps use a single-hue light->dark ramp; the front-minus-steady panel
uses a diverging pair with a neutral midpoint, since its zero is meaningful.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tratt.core import rocking  # noqa: E402
from tratt.cube import front_and_full  # noqa: E402
from tratt.drive import Drive  # noqa: E402
from tratt.observe import channel_stats  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2 = "#0b0b0b", "#52514e"
SEQ = LinearSegmentedColormap.from_list("seq", ["#f2f6fc", "#2a78d6", "#10315c"])
DIV = LinearSegmentedColormap.from_list("div", ["#1b5fa8", "#e8e8e6", "#c04a1f"])

mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.dpi": 200,
})


def compute(alpha, S, alpha_ac, rise=0.0, n_z=1024, n_eta=1501, n_t=401, t_max=1.4, span=2.6):
    d = Drive(S=S, alpha_ac=alpha_ac, rise=rise)
    eta = np.linspace(-span * d.q, span * d.q, n_eta)
    times = np.linspace(0.0, t_max, n_t)
    I_front, I_full = front_and_full(alpha, d, eta, times, n_z=n_z)
    I_static = rocking(alpha, eta, phi=None, n_z=n_z)
    return d, eta, times, I_front, I_full, I_static


def figure(alpha, S, alpha_ac, out="figures", **kw):
    d, eta, times, If, Is, I0 = compute(alpha, S, alpha_ac, **kw)
    x = eta / d.q  # channel units: satellites sit at integers
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)

    # A: the time-resolved rocking curve
    m = ax[0, 0].pcolormesh(x, times, If, cmap=SEQ, shading="auto")
    ax[0, 0].axhline(1.0, color=INK, lw=1.2, ls="--")
    ax[0, 0].text(x.min() * 0.95, 1.02, r"wavefront exit  $\tilde t=1$", va="bottom", fontsize=9)
    fig.colorbar(m, ax=ax[0, 0], label=r"$I_{\rm front}$")
    ax[0, 0].set(xlabel=r"$\eta / q$  (channel index)", ylabel=r"$\tilde t = t / t_{\rm tr}$",
                 title="A  Time-resolved rocking curve")

    # B: departure from the same-phase steady reference
    diff = If - Is
    lim = float(np.max(np.abs(diff)))
    m = ax[0, 1].pcolormesh(x, times, diff, cmap=DIV, shading="auto",
                            norm=TwoSlopeNorm(0.0, -lim, lim))
    ax[0, 1].axhline(1.0, color=INK, lw=1.2, ls="--")
    fig.colorbar(m, ax=ax[0, 1], label=r"$I_{\rm front}-I_{\rm full}$")
    ax[0, 1].set(xlabel=r"$\eta / q$", ylabel=r"$\tilde t$",
                 title="B  Departure from the same-phase steady state")

    # C: channel populations, showing the transient excursion
    colors = {0: BLUE, 1: ORANGE, 2: AQUA}
    for n in (0, 1, 2):
        for sgn, ls in ((1, "-"), (-1, (0, (4, 2)))):
            if n == 0 and sgn == -1:
                continue
            c = channel_stats(eta, If, d, sgn * n)
            ax[1, 0].plot(times, c.p, linestyle=ls, color=colors[n], lw=2,
                          label=f"n={sgn * n:+d}" if n else "n=0")
        cs = channel_stats(eta, Is, d, n)
        ax[1, 0].plot(times, cs.p, linestyle=":", color=colors[n], lw=1.2, alpha=0.8)
        k = int(np.argmax(channel_stats(eta, If, d, n).p))
        ax[1, 0].annotate(f"n={n}", (times[k], channel_stats(eta, If, d, n).p[k]),
                          color=colors[n], fontsize=9, fontweight="bold",
                          xytext=(4, 4), textcoords="offset points")
    ax[1, 0].axvline(1.0, color=INK, lw=1.2, ls="--")
    ax[1, 0].set(xlabel=r"$\tilde t$", ylabel=r"$p_n$  (integrated intensity)",
                 title="C  Channel populations (dotted: steady reference)")
    ax[1, 0].legend(frameon=False, ncol=3, fontsize=9)

    # D: snapshots, so the map in A is readable as curves
    for tt, a_ in zip((0.25, 0.5, 0.75, 1.0), (0.3, 0.5, 0.75, 1.0)):
        k = int(np.argmin(np.abs(times - tt)))
        ax[1, 1].plot(x, If[k], color=BLUE, alpha=a_, lw=1.6, label=rf"$\tilde t$={tt:g}")
    ax[1, 1].plot(x, I0, color=INK2, lw=1.2, ls="--", label="S=0 (no acoustic)")
    ax[1, 1].set(xlabel=r"$\eta / q$", ylabel="I", title="D  Snapshots")
    ax[1, 1].legend(frameon=False, fontsize=9)

    fig.suptitle(rf"Laue,  $\alpha$={alpha:g},  S={S:g},  $\alpha_{{ac}}$={alpha_ac:g},  "
                 rf"rise={kw.get('rise', 0.0):g}", fontsize=12)
    Path(out).mkdir(exist_ok=True)
    p = Path(out) / f"panel_a{alpha:g}_S{S:g}_aac{alpha_ac:g}.png"
    fig.savefig(p)
    plt.close(fig)
    print("wrote", p)
    return p


if __name__ == "__main__":
    figure(1.3, 0.5, 5.0)
    figure(6.0, 1.0, 5.0)
