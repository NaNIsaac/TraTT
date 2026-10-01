"""Three satellites under one front, each at its own stage. Run: python explore/make_ladder_map.py

Laue, S = 3 (J_1 = 0.339, J_2 = 0.486, J_3 = 0.309), alpha = 3, alpha_ac = 60.
Panels a-c: I(eta, t) around n = 1, 2, 3 on a shared colour scale and shared t axis.
Panel d: centre cuts with sin^2(pi alpha |J_n| t). Dashed lines mark the moments
where n = 2 is at its maximum and where it is dark (notebook Entry 24).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

INK, INK2 = "#0b0b0b", "#52514e"
HUES = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a"}
MARK = {1: "o", 2: "s", 3: "^"}
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})


def main(out="figures/poster_ladder_map.png"):
    S, alpha, aac = 3.0, 3.0, 60.0
    d = Drive(S=S, alpha_ac=aac)
    q = d.q
    deta = np.linspace(-40, 40, 321)
    t = np.linspace(0.02, 1.0, 99) + 0.0037
    maps = {n: rocking_cube(alpha, d, n * q + deta, t, n_z=int(60 * aac)) for n in (1, 2, 3)}
    vmax = max(m.max() for m in maps.values())
    t_max2, t_dark2 = 1 / (2 * alpha * jv(2, S)), 1 / (alpha * jv(2, S))

    fig, ax = plt.subplots(1, 4, figsize=(16, 4.6), constrained_layout=True,
                           gridspec_kw=dict(width_ratios=[1, 1, 1, 1.1]))
    k0 = np.argmin(np.abs(deta))
    for i, n in enumerate((1, 2, 3)):
        an = alpha * abs(jv(n, S))
        im = ax[i].imshow(maps[n], origin="lower", aspect="auto", cmap="Blues", vmin=0, vmax=vmax,
                          extent=[deta[0], deta[-1], t[0], t[-1]])
        ax[i].set(title=rf"{'abc'[i]}  $n={n}$:  $\alpha_{n}$={an:.2f}", xlabel=rf"$\eta-{n}q$")
        for tt in (t_max2, t_dark2):
            ax[i].axhline(tt, color=INK2, lw=0.8, ls="--")
    ax[0].set_ylabel(r"$\tilde t=t/t_{\rm tr}$")
    fig.colorbar(im, ax=ax[:3], shrink=0.85, label=r"$I(\eta,\tilde t)$")

    tf = np.linspace(0, 1, 400)
    for n in (1, 2, 3):
        an = alpha * abs(jv(n, S))
        ax[3].plot(tf, np.sin(np.pi * an * tf) ** 2, color=HUES[n], lw=1.6)
        ax[3].plot(t, maps[n][:, k0], MARK[n], ms=3.5, color=HUES[n], mec="white", mew=0.5,
                   ls="none", label=rf"$n={n}$")
    for tt, lab in ((t_max2, "n=2 maximum"), (t_dark2, "n=2 dark")):
        ax[3].axvline(tt, color=INK2, lw=0.8, ls="--")
        ax[3].text(tt + 0.01, 1.08, lab, fontsize=8, color=INK2)
    ax[3].set(title=r"d  Centre cuts; lines: $\sin^2(\pi\alpha|J_n(S)|\tilde t)$", xlabel=r"$\tilde t$",
              ylabel=r"$I(nq,\tilde t)$", ylim=(-0.03, 1.28))
    for s in ("top", "right"):
        ax[3].spines[s].set_visible(False)
    ax[3].legend(fontsize=8, frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle(rf"One front, three satellites at different stages: Laue, $S$={S:g}, $\alpha$={alpha:g},"
                 rf" $\alpha_{{\rm ac}}$={aac:g}", fontsize=10)
    fig.savefig(ROOT / out)
    for tt in (t_max2, t_dark2):
        k = np.argmin(np.abs(t - tt))
        print(f"t={t[k]:.3f}: " + "  ".join(f"n={n}: centre {maps[n][k, k0]:.3f}, "
              f"window {np.trapezoid(maps[n][k], deta):.2f}" for n in (1, 2, 3)))
    print("wrote", out)


if __name__ == "__main__":
    main()
