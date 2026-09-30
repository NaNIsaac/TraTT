"""The two transient features, as one figure. Run: python explore/make_feature_figure.py

Panels are grouped so no axis carries two different measures: the asymmetry
A(t) and the intensity that makes it meaningful are stacked, sharing t, rather
than crammed onto a twin axis. Only three hues are used (slots 1-3), every line
is directly labeled, and the channel sign is carried by line style, so identity
never rests on color alone.

Everything on the right-hand column is computed in the resolved-channel regime
2*pi*alpha < q. Outside it the windows are not satellites and the numbers are
artifacts -- the figure states the criterion rather than leaving it implicit.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tratt.cube import front_and_full, rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402
from tratt.observe import channel_stats  # noqa: E402

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2 = "#0b0b0b", "#52514e"
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})


def channels(alpha, d, eta, times, n_z, ns):
    If, Is = front_and_full(alpha, d, eta, times, n_z=n_z)
    return ({n: channel_stats(eta, If, d, n) for n in ns},
            {n: channel_stats(eta, Is, d, n) for n in ns})


def main(out="figures/feature_overshoot_asymmetry.png"):
    fig, ax = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)

    # ---- column 1: asymmetry of the reference case, with its validity weight
    d = Drive(S=0.5, alpha_ac=5.0)
    eta = np.linspace(-3 * d.q, 3 * d.q, 1501)
    t = np.linspace(0.0, 1.6, 321)
    cf, cs = channels(1.3, d, eta, t, 768, (1, -1))
    Af = (cf[1].p - cf[-1].p) / (cf[1].p + cf[-1].p)
    As = (cs[1].p - cs[-1].p) / (cs[1].p + cs[-1].p)
    ax[0, 0].plot(t, Af, color=ORANGE, lw=2)
    ax[0, 0].plot(t, As, color=INK2, lw=1.4, linestyle=":")
    ax[0, 0].annotate("traveling front", (0.35, Af[t.searchsorted(0.35)]), color=ORANGE,
                      fontweight="bold", fontsize=9, xytext=(4, 6), textcoords="offset points")
    ax[0, 0].annotate("steady reference", (1.3, As[t.searchsorted(1.3)]), color=INK2,
                      fontsize=9, xytext=(-30, -16), textcoords="offset points")
    ax[0, 0].axvline(1.0, color=INK, lw=1.2, linestyle="--")
    ax[0, 0].axhline(0.0, color=INK2, lw=0.8)
    ax[0, 0].set(ylabel=r"$A_1=(p_{+1}-p_{-1})/(p_{+1}+p_{-1})$",
                 title=r"a  Traveling front breaks $\pm n$ symmetry")

    ax[1, 0].plot(t, cf[1].p + cf[-1].p, color=ORANGE, lw=2)
    ax[1, 0].axvline(1.0, color=INK, lw=1.2, linestyle="--")
    ax[1, 0].set(xlabel=r"$\tilde t=t/t_{\rm tr}$", ylabel=r"$p_{+1}+p_{-1}$",
                 title="b  Intensity carrying that asymmetry")
    ax[1, 0].annotate("early peak of |A| sits\nwhere the signal is weak",
                      (0.15, (cf[1].p + cf[-1].p)[t.searchsorted(0.15)]), fontsize=8,
                      color=INK2, xytext=(10, 30), textcoords="offset points",
                      arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))

    # ---- column 2: how both features depend on alpha
    alphas = np.array([0.2, 0.5, 1.0, 2.0, 3.0])
    maxA, rmsA, steadyA = [], [], []
    for a in alphas:
        cf2, cs2 = channels(a, d, eta, t, 768, (1, -1))
        A = (cf2[1].p - cf2[-1].p) / (cf2[1].p + cf2[-1].p)
        S_ = (cs2[1].p - cs2[-1].p) / (cs2[1].p + cs2[-1].p)
        pre = t < 1.0
        maxA.append(np.max(np.abs(A[pre])))
        rmsA.append(np.sqrt(np.mean(A[pre] ** 2)))
        steadyA.append(np.max(np.abs(S_)))
    ax[0, 1].plot(alphas, maxA, "o-", color=ORANGE, lw=2, ms=7)
    ax[0, 1].plot(alphas, rmsA, "s--", color=ORANGE, lw=1.6, ms=6, alpha=0.7)
    ax[0, 1].plot(alphas, steadyA, "^-", color=INK2, lw=1.6, ms=6)
    for y, lab, col in ((maxA, r"max$|A|$, $\tilde t<1$", ORANGE),
                        (rmsA, r"rms$|A|$, $\tilde t<1$", ORANGE),
                        (steadyA, r"max$|A|$, steady", INK2)):
        ax[0, 1].annotate(lab, (alphas[-1], y[-1]), color=col, fontsize=9,
                          xytext=(-4, 6), textcoords="offset points", ha="right")
    ax[0, 1].set(xlabel=r"$\alpha=L/\Lambda_{\rm ext}$", ylabel=r"$|A_1|$", ylim=(0, None),
                 title=r"c  Asymmetry is nearly $\alpha$-independent")

    # ---- overshoot in the resolved regime
    d2 = Drive(S=1.0, alpha_ac=20.0)
    eta2 = np.linspace(-2.5 * d2.q, 2.5 * d2.q, 2001)
    t2 = np.linspace(0.0, 1.2, 241)
    ratios = {n: [] for n in (0, 1, 2)}
    al2 = np.array([1.0, 2.0, 4.0, 6.0])
    for a in al2:
        cf3, cs3 = channels(a, d2, eta2, t2, 1024, (0, 1, 2))
        for n in (0, 1, 2):
            ratios[n].append(cf3[n].p.max() / np.mean(cs3[n].p))
    for n, col in zip((0, 1, 2), (BLUE, ORANGE, AQUA)):
        ax[1, 1].plot(al2, ratios[n], "o-", color=col, lw=2, ms=7, label=f"n={n}")
        ax[1, 1].annotate(f"n={n}", (al2[-1], ratios[n][-1]), color=col, fontweight="bold",
                          fontsize=9, xytext=(4, 0), textcoords="offset points")
    ax[1, 1].axhline(1.0, color=INK2, lw=1.0, linestyle=":")
    ax[1, 1].set(xlabel=r"$\alpha$", ylabel=r"$\max_t p_n \,/\, \langle p_n^{\rm full}\rangle$",
                 title=r"d  Transient overshoot ($\alpha_{ac}$=20, channels resolved)")
    ax[1, 1].legend(frameon=False, fontsize=9, ncol=3)

    # ---- column 3: the overshoot itself, and its independence of the rise time
    cf4, cs4 = channels(2.0, d2, eta2, t2, 1024, (0, 1, 2))
    for n, col in zip((0, 1, 2), (BLUE, ORANGE, AQUA)):
        norm = np.mean(cs4[n].p)
        ax[0, 2].plot(t2, cf4[n].p / norm, color=col, lw=2)
        k = int(np.argmax(cf4[n].p))
        ax[0, 2].annotate(f"n={n}", (t2[k], cf4[n].p[k] / norm), color=col, fontweight="bold",
                          fontsize=9, xytext=(4, 2), textcoords="offset points")
    ax[0, 2].axhline(1.0, color=INK2, lw=1.0, linestyle=":")
    ax[0, 2].axvline(1.0, color=INK, lw=1.2, linestyle="--")
    ax[0, 2].set(xlabel=r"$\tilde t$", ylabel=r"$p_n/\langle p_n^{\rm full}\rangle$",
                 title=r"e  Overshoot, $\alpha$=2 (dotted: steady level)")

    for rise, alpha_v in zip((0.0, 0.3, 0.6), (1.0, 0.65, 0.4)):
        dr = Drive(S=1.0, alpha_ac=20.0, rise=rise)
        tr = np.linspace(0.0, 1.3 + rise, 241)
        I = rocking_cube(2.0, dr, eta2, tr, n_z=1024)
        Isr = rocking_cube(2.0, dr, eta2, tr, n_z=1024, steady=True)
        p = channel_stats(eta2, I, dr, 1).p
        ax[1, 2].plot(tr, p / np.mean(channel_stats(eta2, Isr, dr, 1).p),
                      color=ORANGE, lw=2, alpha=alpha_v, label=f"rise={rise:g}")
    ax[1, 2].axhline(1.0, color=INK2, lw=1.0, linestyle=":")
    ax[1, 2].set(xlabel=r"$\tilde t$", ylabel=r"$p_{+1}/\langle p_{+1}^{\rm full}\rangle$",
                 title="f  Overshoot survives a finite transducer rise")
    ax[1, 2].legend(frameon=False, fontsize=9)

    Path(out).parent.mkdir(exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
