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
sys.path.insert(0, str(ROOT / "src"))

from tratt.core import required_nz, rocking  # noqa: E402
from tratt.drive import Drive  # noqa: E402
OUT = ROOT / "figures" / "poster"
for _f in sorted((ROOT / "poster" / "fonts").glob("*.ttf")):     # same typeface as the poster text
    mpl.font_manager.fontManager.addfont(str(_f))
FONT = "IBM Plex Sans"

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
        # Text in the poster's own face; math in STIX, the same pairing poster/build.py uses, so
        # an axis label and a panel equation are set identically. An all-Plex math fontset has no
        # \cdot and no \sim and falls back to DejaVu for them, and its \int is a text-sized glyph.
        "font.family": FONT, "mathtext.fontset": "stix",
    })


def finish(fig, ax, name):
    for a in np.atleast_1d(ax):
        a.grid(color=GRID, lw=1.0)
        a.set_axisbelow(True)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png")
    print("wrote", OUT / f"{name}.png")


def f0():
    """Schematic for half a wide card (text sits on the other half), so FIG_W / 2 wide.

    Everything here is computed, like the other four figures: the profile inside the slab is
    `Drive.phi(t)` sampled in depth, and the curve on the right is a real `rocking()` of that
    same snapshot, turned on its side so each satellite peak faces the beam it belongs to.
    The earlier version drew sine curves wiggling along the transverse direction, i.e. the
    modulation varying across the beam rather than with depth -- the opposite of the drive.
    """
    style()
    alpha, S, alpha_ac, t = 0.6, 1.5, 5.0, 0.62
    d = Drive(S=S, alpha_ac=alpha_ac)

    fig, (ax, axr) = plt.subplots(
        1, 2, figsize=(FIG_W / 2, 3.2), constrained_layout=True,
        gridspec_kw=dict(width_ratios=[2.5, 1.0], wspace=0.0))

    # ---- left: the crystal in cross-section, with the real phase profile inside it
    ax.set(xlim=(-3.3, 9.6), ylim=(-0.3, 9.6))
    ax.axis("off")
    x0, x1, y0, y1 = 0.0, 6.0, 2.1, 8.1            # slab; depth zeta runs left to right
    xf = x0 + t * (x1 - x0)                         # acoustic front
    ymid = 0.5 * (y0 + y1)

    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fc="#f4f7fb", ec=INK2, lw=1.6, zorder=1))
    ax.add_patch(plt.Rectangle((x0, y0), xf - x0, y1 - y0, fc="#e7effa", ec="none", zorder=1))

    zz = np.linspace(0.0, 1.0, 1200)
    phi = d.phi(t)(zz)
    amp = 0.33 * (y1 - y0) / max(abs(phi).max(), 1e-9)
    ax.plot(x0 + zz * (x1 - x0), ymid + amp * phi, color=BLUE, lw=2.6, zorder=3,
            solid_capstyle="round")
    ax.plot([xf, x1], [ymid, ymid], color=BLUE, lw=2.6, zorder=3, solid_capstyle="round")
    ax.plot([xf, xf], [y0, y1], color=BLUE, lw=2.4, ls=(0, (4, 3)), zorder=4)

    # transducer on the lower part of the entrance face, so the beam enters clear of it
    ty = y0 + 0.58 * (y1 - y0)
    ax.add_patch(plt.Rectangle((x0 - 0.44, y0), 0.44, ty - y0, fc=INK2, ec="none", zorder=2))
    ax.text(x0 - 0.62, 0.5 * (y0 + ty), "transducer", ha="right", va="center",
            fontsize=SMALL, color=INK2)

    # front marker
    ax.annotate("", xy=(xf + 1.25, y1 + 0.70), xytext=(xf, y1 + 0.70),
                arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.2, mutation_scale=22))
    ax.text(xf + 1.45, y1 + 0.70, r"$v_s$", color=BLUE, va="center", fontsize=BIG)

    # one dimension row: the slab thickness, with the filled part picked out on the same line
    yd = y0 - 0.72
    ax.annotate("", xy=(x1, yd), xytext=(x0, yd),
                arrowprops=dict(arrowstyle="<|-|>", color=INK2, lw=1.2, mutation_scale=13))
    for x in (x0, xf, x1):
        ax.plot([x, x], [yd - 0.22, yd + 0.22], color=INK2, lw=1.1)
    ax.plot([x0, xf], [yd, yd], color=INK, lw=3.2, solid_capstyle="butt")
    ax.text(x1 + 0.22, yd, r"$L$", ha="left", va="center", fontsize=SMALL, color=INK2)
    ax.text(0.5 * (x0 + xf), yd - 0.38, r"filled depth $v_s t$", ha="center", va="top",
            fontsize=SMALL, color=INK)

    # beams: in above the transducer, forward and diffracted out
    ax.annotate("", xy=(x0 + 0.08, ymid + 1.25), xytext=(-2.8, ymid + 2.85),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=2.8, mutation_scale=26), zorder=5)
    ax.text(-2.9, ymid + 3.2, "X-ray", fontsize=BIG, color=INK)
    ax.annotate("", xy=(9.0, ymid + 1.65), xytext=(x1, ymid + 0.1),
                arrowprops=dict(arrowstyle="-|>", color=INK2, lw=2.4, mutation_scale=22))
    ax.text(9.2, ymid + 1.8, r"$D_0$", va="center", fontsize=BIG, color=INK2)
    ax.annotate("", xy=(9.0, ymid - 1.65), xytext=(x1, ymid - 0.1),
                arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=3.0, mutation_scale=22))
    ax.text(9.2, ymid - 1.8, r"$D_h$", va="center", fontsize=BIG, color=ORANGE)

    # ---- right: the same snapshot, actually solved, eta vertical so peaks face the beam
    eta = np.linspace(-2.6 * d.q, 2.6 * d.q, 4001)
    I = rocking(alpha, eta, phi=d.phi(t), n_z=required_nz(alpha, alpha_ac))
    axr.plot(I, eta / d.q, color=ORANGE, lw=2.2)
    axr.fill_betweenx(eta / d.q, 0, I, color=ORANGE, alpha=0.16)
    axr.set(xlim=(0, 1.55 * I.max()), ylim=(-2.75, 2.75), xticks=[], yticks=[-2, -1, 0, 1, 2])
    axr.set_xlabel(r"$I(\eta)$", fontsize=SMALL, labelpad=1)
    axr.set_ylabel(r"$\eta / q$", fontsize=SMALL, labelpad=0)
    axr.tick_params(labelsize=SMALL, length=4)
    for s in ("top", "right", "bottom"):
        axr.spines[s].set_visible(False)
    axr.grid(False)
    for n in (0, 1, 2):
        axr.annotate(f"$n={n}$", (float(np.interp(n * d.q, eta, I)), n), fontsize=SMALL,
                     color=ORANGE, va="center", xytext=(6, 0), textcoords="offset points")

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "f0_scheme.png", dpi=400)  # half the card width, so twice the dpi to print like the rest
    print("wrote", OUT / "f0_scheme.png")


def f1():
    import make_poster_figure as M
    style()
    fig, ax = plt.subplots(1, 2, figsize=(FIG_W, 5.9), constrained_layout=True, sharey=True)
    M.collapse(ax[0], "laue", lambda s: np.sin(np.pi * s) ** 2, r"$\sin^2(\pi\alpha_n\tilde t)$",
               r"Laue: $\sin^2(\pi\alpha_n\tilde t)$", small=SMALL - 4, ms=9)
    M.collapse(ax[1], "bragg", lambda s: np.tanh(np.pi * s) ** 2, r"$\tanh^2(\pi\alpha_n\tilde t)$",
               r"Bragg: $\tanh^2(\pi\alpha_n\tilde t)$", small=SMALL - 4, ms=9)
    ax[0].get_legend().remove()                      # one legend for both: the cases are the same
    h, lab = ax[1].get_legend_handles_labels()
    ax[1].legend(h, lab[:-1] + ["closed form (line)"], frameon=False, loc="lower right", fontsize=SMALL - 2)
    for a in ax:
        a.set_xlabel(r"$\alpha_n\tilde t$,   $\alpha_n=\alpha\left|J_n(S)\right|$")
    ax[1].set_ylabel("")
    ax[1].tick_params(labelleft=False)          # shared y axis: one set of tick labels
    finish(fig, ax, "f1_collapse")


def f2():
    import make_map_figure as M
    style()
    r = M.compute()
    deta, t, I, an = r["deta"], r["t"], r["I"], r["an"]
    fig, ax = plt.subplots(1, 2, figsize=(FIG_W, 5.9), constrained_layout=True,
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
    fig, ax = plt.subplots(1, 2, figsize=(FIG_W, 5.9), constrained_layout=True)
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
    fig, ax = plt.subplots(figsize=(FIG_W, 5.0), constrained_layout=True)
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
