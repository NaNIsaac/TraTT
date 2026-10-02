"""Poster-only figures, drawn at print size for one A0 column. Run: python explore/make_poster_figs.py [f0 f1 f2 f3 f4]

Width FIG_W is the inner width of a poster card (poster/build.py SPEC), so a font size here is the
printed size. Panels reuse the notebook figure code; nothing is recomputed differently.

f0  schematic: Laue slab, acoustic front filling the crystal, incident / transmitted / satellite beams
f1  filled-depth mapping: Laue sin^2 and Bragg tanh^2 collapse      (make_poster_figure.collapse; Entries 1, 9)
f2  time-resolved rocking curve near n = 1 and its centre cut        (make_map_figure.compute; Entries 22, 23)
f3  path memory: transient-only satellite, up vs down a pulse       (make_memory_figure.panel_a/b; Entries 28, 30)
f4  population error of the independent-channel picture vs R        (make_phase_figure.merged; Entry 54)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "explore"))
OUT = ROOT / "figures" / "poster"

FIG_W = 14.4          # inches = poster card inner width (~366 mm)
BIG, SMALL = 24, 20   # printed pt
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"


def style():
    mpl.rcParams.update({
        "font.size": BIG, "axes.titlesize": BIG, "axes.labelsize": BIG, "xtick.labelsize": SMALL,
        "ytick.labelsize": SMALL, "legend.fontsize": SMALL, "lines.linewidth": 3.2, "lines.markersize": 10,
        "axes.linewidth": 1.4, "xtick.major.width": 1.4, "ytick.major.width": 1.4,
        "xtick.major.size": 7, "ytick.major.size": 7, "axes.titlepad": 12,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2,
        "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.dpi": 200,
        "font.family": "Arial", "mathtext.fontset": "dejavusans",
    })


def finish(fig, ax, name):
    for a in np.atleast_1d(ax):
        a.grid(color=GRID, lw=1.0)
        a.set_axisbelow(True)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png")
    print("wrote", OUT / f"{name}.png")


def f0():
    """Schematic for half a wide card (text sits on the other half), so FIG_W / 2 wide."""
    style()
    fig, ax = plt.subplots(figsize=(FIG_W / 2, 3.6), constrained_layout=True)
    ax.set(xlim=(0, 14), ylim=(0, 9.2))
    ax.axis("off")
    x0, x1, y0, y1 = 3.2, 7.2, 1.6, 7.6                       # crystal slab (Laue: beams cross it)
    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fc="#eef2f8", ec=INK2, lw=1.6))
    xf = x0 + 0.62 * (x1 - x0)                                  # acoustic front inside the crystal
    yy = np.linspace(y0, y1, 400)
    for k in range(6):
        xs = x0 + 0.05 + k * (xf - x0 - 0.1) / 5
        ax.plot(xs + 0.12 * np.sin(2 * np.pi * (yy - y0) / 1.5), yy, color=BLUE, lw=1.6, alpha=0.8)
    ax.plot([xf, xf], [y0, y1], color=BLUE, lw=3, ls="--")
    ax.annotate("", xy=(xf + 1.6, 8.3), xytext=(xf - 0.2, 8.3),
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.5, mutation_scale=26))
    ax.text(xf + 1.8, 8.3, r"front, $v_s$", color=BLUE, va="center", fontsize=BIG)
    ax.add_patch(plt.Rectangle((x0 - 0.35, y0), 0.35, y1 - y0, fc=INK2, ec="none"))
    ax.text(x0 - 0.2, 0.75, "transducer", ha="center", fontsize=SMALL, color=INK2)
    ax.annotate("", xy=(x0 + 0.1, 4.6), xytext=(0.2, 6.8),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=3, mutation_scale=28))
    ax.text(0.1, 7.3, "X-ray", fontsize=BIG)
    ax.annotate("", xy=(12.4, 7.4), xytext=(x1, 4.6),       # forward beam
                arrowprops=dict(arrowstyle="-|>", color=INK2, lw=2.6, mutation_scale=24))
    ax.text(12.6, 7.4, r"$D_0$", va="center", fontsize=BIG, color=INK2)
    for n in (-2, -1, 0, 1, 2):                              # diffracted beam D_h and its satellites
        y = 3.0 - 0.75 * n
        ax.annotate("", xy=(12.4, y), xytext=(x1, 4.6),
                    arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=3.2 if n == 0 else 1.8,
                                    mutation_scale=24 if n == 0 else 16, alpha=1 if n == 0 else 0.75))
    ax.text(12.6, 3.0, r"$D_h$", va="center", fontsize=BIG, color=ORANGE)
    ax.text(13.9, 0.75, r"satellites $n=\pm1,\pm2$", ha="right", fontsize=SMALL, color=ORANGE)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "f0_scheme.png")
    print("wrote", OUT / "f0_scheme.png")


def f1():
    import make_poster_figure as M
    style()
    fig, ax = plt.subplots(1, 2, figsize=(FIG_W, 5.4), constrained_layout=True)
    M.collapse(ax[0], "laue", lambda s: np.sin(np.pi * s) ** 2, r"$\sin^2(\pi\alpha_n\tilde t)$",
               r"Laue: $\sin^2(\pi\alpha_n\tilde t)$", small=SMALL - 4, ms=9)
    M.collapse(ax[1], "bragg", lambda s: np.tanh(np.pi * s) ** 2, r"$\tanh^2(\pi\alpha_n\tilde t)$",
               r"Bragg: $\tanh^2(\pi\alpha_n\tilde t)$", small=SMALL - 4, ms=9)
    ax[0].get_legend().remove()                      # one legend for both: the cases are the same
    h, lab = ax[1].get_legend_handles_labels()
    ax[1].legend(h, lab[:-1] + ["closed form (line)"], frameon=False, loc="lower right", fontsize=SMALL - 2)
    for a in ax:
        a.set_xlabel(r"$\alpha_n\tilde t$,   $\alpha_n=\alpha|J_n(S)|$")
    ax[1].set_ylabel("")
    finish(fig, ax, "f1_collapse")


def f2():
    import make_map_figure as M
    style()
    r = M.compute()
    deta, t, I, an = r["deta"], r["t"], r["I"], r["an"]
    fig, ax = plt.subplots(1, 2, figsize=(FIG_W, 5.4), constrained_layout=True,
                           gridspec_kw=dict(width_ratios=[1.25, 1]))
    im = ax[0].imshow(I, origin="lower", aspect="auto", cmap="Blues", vmin=0, vmax=I.max(),
                      extent=[deta[0], deta[-1], t[0], t[-1]])
    ax[0].set(xlabel=r"$\eta-q$", ylabel=r"$\tilde t=t/t_{\rm tr}$", title=r"$I(\eta,\tilde t)$, full TT")
    fig.colorbar(im, ax=ax[0], shrink=0.9, pad=0.02)
    k0 = np.argmin(np.abs(deta))
    tf = np.linspace(0, 1, 400)
    ax[1].plot(tf, np.sin(np.pi * an * tf) ** 2, color=INK, lw=2.4, label=r"$\sin^2(\pi\alpha_1\tilde t)$")
    ax[1].plot(t, I[:, k0], "o", color=BLUE, mec="white", mew=0.8, label="full TT")
    ax[1].set(xlabel=r"$\tilde t$", ylabel=r"$I(\eta=q,\tilde t)$", ylim=(-0.03, 1.1),
              title=r"centre: TT vs $\sin^2(\pi\alpha_1\tilde t)$")
    for s in ("top", "right"):
        ax[0].spines[s].set_visible(True)
    finish(fig, ax[1], "f2_fringes")


def f3():
    import make_memory_figure as M
    style()
    fig, ax = plt.subplots(1, 2, figsize=(FIG_W, 5.4), constrained_layout=True)
    M.panel_a(ax[0], small=SMALL - 2, title=r"transient-only: $J_1(S_f)=0$")
    M.panel_b(ax[1], small=SMALL - 4, title="same strain, up vs down")
    for line in ax[0].lines + ax[1].lines:
        if line.get_marker() not in ("None", None, ""):
            line.set_markersize(8)
    ax[0].texts[0].set(position=(2.93, 0.66), ha="right", va="top")
    ax[0].legend(frameon=False, loc="upper left", fontsize=SMALL)
    ax[1].set_ylim(-0.02, 0.76)
    ax[1].legend(frameon=False, loc="upper center", ncol=2, fontsize=SMALL - 2, columnspacing=1.0, handletextpad=0.3)
    # which branch is which: every loop starts at the origin, so the lower branch is the rise
    ax[1].text(2.05, 0.09, "up", color=BLUE, fontsize=SMALL, fontweight="bold")
    ax[1].text(0.02, 0.475, "down", color=BLUE, fontsize=SMALL, fontweight="bold")
    finish(fig, ax, "f3_memory")


def f4():
    import make_phase_figure as M
    style()
    P = M.merged()
    fig, ax = plt.subplots(figsize=(FIG_W, 4.5), constrained_layout=True)
    Ss = sorted(set(P[:, 2]))
    hues = plt.cm.viridis(np.linspace(0, 0.9, len(Ss)))
    for c, S in zip(hues, Ss):
        m = (P[:, 2] == S) & (P[:, 1] >= 2)
        ax.plot(P[m, 0] / P[m, 1], P[m, 3], "o", ms=8, color=c, mec="white", mew=0.5, ls="none", label=f"S={S:g}")
    for e, lab in ((0.05, "5 %"), (0.2, "20 %"), (0.5, "50 %")):
        ax.axhline(e, color=INK2, lw=1.4, ls=":")
        ax.text(1.3e-3, e * 1.12, lab, fontsize=SMALL, color=INK2)
    ax.set(xscale="log", yscale="log", ylim=(5e-5, 1.5), xlim=(1e-3, 6),
           xlabel=r"$R=\Lambda_{ac}/\Lambda_P$", ylabel="population error $E_{pop}$")
    # nine S values: a discrete colour bar reads better than a nine-entry legend
    cmap = mpl.colors.ListedColormap(hues)
    norm = mpl.colors.BoundaryNorm(np.arange(len(Ss) + 1) - 0.5, len(Ss))
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.01, aspect=25)
    cb.set_ticks(range(len(Ss)), labels=[f"{S:g}" for S in Ss])
    cb.set_label("modulation index $S$")
    cb.outline.set_visible(False)
    finish(fig, ax, "f4_regime")


if __name__ == "__main__":
    for name in (sys.argv[1:] or ["f0", "f1", "f2", "f3", "f4"]):
        globals()[name]()
