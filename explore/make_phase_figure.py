"""Regime map of the resolved Floquet-channel picture. Run: python explore/make_phase_figure.py

Uniform travelling wave, filled crystal, Laue. E_pop = order-window population error of the
independent-channel picture (Entry 43). Phase-averaged: exact sum_m |b_m|^2 where computed
(cache_e43_exact, Entry 54), else 16 phases (cache_e43_ph16, Entry 53), else 4 phases (cache_e43_ph4).
a  E_pop over (alpha, alpha_ac) at S = 1.5, with E_pop = 0.05 / 0.2 / 0.5 contours.
b  E_pop vs R for all S (alpha_ac >= 2), phase-averaged (filled) and single phase t = 0 (open);
   thresholds 0.05 / 0.2 / 0.5. Boundaries in R depend on S (Entry 54 table).
c  total diffracted intensity / independent prediction vs R (16 phases where available, else 4).
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
mpl.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": "white",
                     "savefig.facecolor": "white", "savefig.dpi": 200})


def load(name):
    rows = [json.loads(open(f).read()) for f in glob.glob(str(ROOT / "explore" / "nb" / name / "*.json"))]
    return np.array([[r["alpha"], r["aac"], r["S"], r["E_pop"], r["tot"] / max(r["tot_ind"], 1e-30)]
                     for r in rows if r["S"] > 0])


def merged():
    """Best available phase average per point: exact > 16 phases > 4 phases (E_pop); tot ratio 16 > 4."""
    best = {tuple(r[:3]): r.copy() for r in load("cache_e43_ph4")}
    for r in load("cache_e43_ph16"):
        best[tuple(r[:3])] = r.copy()
    for f in glob.glob(str(ROOT / "explore" / "nb" / "cache_e43_exact" / "*.json")):
        r = json.loads(open(f).read())
        k = (r["alpha"], r["aac"], r["S"])
        if k in best:
            best[k][3] = r["E_pop"]
    return np.array(list(best.values()))


def main(out="figures/regime_phase_diagram.png"):
    P, T0 = merged(), load("cache_e43")
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
    S0 = 1.5
    sel = P[:, 2] == S0
    als, aacs = sorted(set(P[sel, 0])), sorted(set(P[sel, 1]))
    M = np.full((len(als), len(aacs)), np.nan)
    for a, q, s, e, _ in P[sel]:
        M[als.index(a), aacs.index(q)] = e
    im = ax[0].pcolormesh(aacs, als, np.log10(M), cmap="Blues", vmin=-3, vmax=0, shading="nearest")
    cs = ax[0].contour(aacs, als, M, levels=[0.05, 0.2, 0.5], colors=INK, linewidths=1.0,
                       linestyles=[":", "--", "-"])
    ax[0].clabel(cs, fmt={0.05: "0.05", 0.2: "0.2", 0.5: "0.5"}, fontsize=8)
    ax[0].set(xscale="log", yscale="log", xlim=(0.45, 45), ylim=(0.045, 11), xlabel=r"$\alpha_{ac}=L/\Lambda_{ac}$",
              ylabel=r"$\alpha=L/\Lambda_P$", title=rf"a  $\log_{{10}}E_{{pop}}$, S={S0} (phase-averaged)")
    fig.colorbar(im, ax=ax[0], shrink=0.8)
    hues = plt.cm.viridis(np.linspace(0, 0.9, 9))
    for c, S in zip(hues, sorted(set(P[:, 2]))):
        for D, mfc in ((P, None), (T0, "white")):
            m = (D[:, 2] == S) & (D[:, 1] >= 2)
            ax[1].plot(D[m, 0] / D[m, 1], D[m, 3], "o", ms=3.5, color=c, mfc=mfc or c, ls="none",
                       label=f"S={S:g}" if mfc is None else None)
        m = (P[:, 2] == S) & (P[:, 1] >= 2)
        ax[2].plot(P[m, 0] / P[m, 1], P[m, 4], "o", ms=3.5, color=c, ls="none", label=f"S={S:g}")
    for e, lab in ((0.05, "A|B"), (0.2, "B|C"), (0.5, "C|D")):
        ax[1].axhline(e, color=INK2, lw=0.8, ls=":")
        ax[1].text(3.5e-3, e * 1.08, lab, fontsize=8, color=INK2)
    ax[1].set(xscale="log", yscale="log", ylim=(5e-5, 2), xlabel=r"$R=\alpha/\alpha_{ac}=\Lambda_{ac}/\Lambda_P$",
              ylabel=r"$E_{pop}$", title=r"b  Breakdown vs R (filled: phase-avg, open: one phase)")
    ax[1].legend(fontsize=7, frameon=False, ncol=3, loc="lower right")
    ax[2].set(xscale="log", ylim=(0, 1.25), xlabel="R", ylabel="total diffracted / independent prediction",
              title="c  Overlapping orders share the forward beam")
    for a in ax:
        a.grid(color=GRID, lw=0.6)
        a.set_axisbelow(True)
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
