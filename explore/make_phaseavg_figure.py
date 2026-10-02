"""Phase-average convergence and the observability of the Bessel-zero residual. Run: python explore/make_phaseavg_figure.py

phaseavg_convergence.png : |<I>_N - <I>_exact| (absolute) and relative, vs number N of uniform acoustic
                           phases, at eta = q, alpha = 1, S = j11, alpha_ac = 20, 20.25, 20.5
                           (exact average = sum_m |b_m(1)|^2, notebook Entry 53).
bessel_residual_eps.png  : alpha_ac = 20 + eps: target population |b_-1|^2, exact phase-averaged exit
                           intensity, and the rest (sum over m != -1), at eta = q; target fraction.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "explore" / "nb"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(NB))
from e53_phaseavg import J11, components, phase_avg_N  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
HUES = ["#2a78d6", "#eb6834", "#1baf7a"]
mpl.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": "white",
                     "savefig.facecolor": "white", "savefig.dpi": 200})


def convergence():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    Ns = np.arange(1, 65)
    for aac, c in zip((20.0, 20.25, 20.5), HUES):
        q = 2 * np.pi * aac
        b, ms, _ = components(1.0, aac, J11, np.array([q]))
        ex = float(np.sum(np.abs(b[0]) ** 2))
        av = np.array([phase_avg_N(b, ms, q, N)[0] for N in Ns])
        ax[0].plot(Ns, np.maximum(np.abs(av - ex), 1e-18), "o-", ms=3, color=c, label=rf"$\alpha_{{ac}}$={aac:g}")
        ax[1].plot(Ns, np.maximum(np.abs(av - ex) / ex, 1e-14), "o-", ms=3, color=c, label=rf"$\alpha_{{ac}}$={aac:g}")
    for a in ax:
        a.axvline(4, color=INK2, lw=0.8, ls=":")
        a.text(4.3, 1e-3 if a is ax[1] else 1e-6, "N = 4 (Entry 43)", fontsize=8, color=INK2)
        a.set(xscale="log", yscale="log", xlabel="number of uniform acoustic phases N")
    ax[0].set(ylabel=r"$|\langle I\rangle_N-\langle I\rangle_{\rm exact}|$ at $\eta=q$",
              title=r"a  Absolute error ($\alpha$=1, $S=j_{1,1}$)")
    ax[1].set(ylabel="relative error", title="b  Relative error")
    ax[0].legend(fontsize=8, frameon=False)
    for a in ax:
        a.grid(color=GRID, lw=0.6)
    fig.savefig(ROOT / "figures" / "phaseavg_convergence.png")


def eps_scan():
    rows = json.load(open(NB / "e53_eps.json"))
    e = np.array([r["eps"] for r in rows])
    tgt, tot = np.array([r["tgt_nom"] for r in rows]), np.array([r["tot_nom"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3), constrained_layout=True)
    for a, lim in ((ax[0], 0.5), (ax[1], 0.01)):
        m = np.abs(e) <= lim
        a.plot(e[m], tot[m], color=INK, lw=1.4, label="exact phase-averaged exit intensity")
        a.plot(e[m], tgt[m], color=HUES[0], lw=1.4, label=r"target $|b_{-1}|^2$ (third-order residual)")
        a.plot(e[m], np.maximum(tot[m] - tgt[m], 1e-16), color=HUES[1], lw=1.2, ls="--",
               label=r"rest $\sum_{m\neq-1}|b_m|^2$")
        a.set(yscale="log", xlabel=r"$\epsilon$ ($\alpha_{ac}=20+\epsilon$)", ylabel=r"intensity at $\eta=q$")
        a.grid(color=GRID, lw=0.6)
    ax[0].set_title(r"a  $\alpha$=1, $S=j_{1,1}$: the residual dominates only at integer $\alpha_{ac}$")
    ax[0].legend(fontsize=8, frameon=False, loc="lower right")
    f = tgt / tot
    axf = ax[1].twinx()  # target fraction (same eta, different quantity) - kept on its own axis, labelled
    axf.plot(e[np.abs(e) <= 0.01], f[np.abs(e) <= 0.01], color=HUES[2], lw=1.2)
    for thr in (0.5, 0.9):
        axf.axhline(thr, color=HUES[2], lw=0.6, ls=":")
    axf.set_ylabel("target fraction of exit intensity", color=HUES[2])
    axf.spines["right"].set_visible(True)
    ax[1].set_title(r"b  zoom: target $\geq$ 50 % for $|\epsilon|\leq0.004$, $\geq$ 90 % for $|\epsilon|\leq0.0014$")
    fig.savefig(ROOT / "figures" / "bessel_residual_eps.png")


if __name__ == "__main__":
    convergence()
    eps_scan()
    print("wrote figures")
