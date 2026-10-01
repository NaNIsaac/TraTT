"""Summary of the coarse anomaly search (notebook Entry 40). Run: python explore/make_search_figure.py

a  integrated satellite relative to the uniform-q control, full TT vs kinematic, for every
   case, against alpha: the kinematic ratio stays at 1 for phase-only profiles (Parseval);
   the dynamical ratio grows with alpha.
b  line width vs the uniform control: nothing goes below the dynamical-narrowing floor
   0.71 x 5.566 (Entry 4).
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
FAM = {"chirp": ("#2a78d6", "o"), "qwall": ("#eb6834", "s"), "pwall": ("#1baf7a", "^"), "ramp": ("#eda100", "D")}
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})


def main(out="figures/search_summary.png"):
    rows = json.load(open(ROOT / "explore" / "nb" / "e40_search.json"))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    rng = np.random.default_rng(0)
    for fam, (c, m) in FAM.items():
        r = [x for x in rows if x["fam"] == fam]
        a = np.array([x["alpha"] for x in r]) * np.exp(0.06 * rng.standard_normal(len(r)))
        ax[0].plot(a, [x["int_ratio"] for x in r], m, ms=5, color=c, ls="none", label=f"{fam} (full TT)")
        ax[0].plot(a, [np.exp(x["kin_dint"]) if x["fam"] != "ramp" else np.nan for x in r], m, ms=5,
                   mfc="white", color=c, ls="none")
        ax[1].plot([x["width_c"] for x in r], [x["width"] for x in r], m, ms=5, color=c, ls="none", label=fam)
    ax[0].axhline(1, color=INK2, lw=0.8, ls=":")
    ax[0].set(xscale="log", yscale="log", xlabel=r"$\alpha$ (jittered)", ylabel="integrated satellite / uniform-q",
              title="a  New spectral weight is dynamical (open: kinematic = 1)")
    ax[0].legend(fontsize=8, frameon=False, loc="upper left")
    w = np.array([3, 300])
    ax[1].plot(w, w, color=INK2, lw=0.8, ls=":")
    ax[1].axhline(0.71 * 5.566, color=INK, lw=1.0, ls="--")
    ax[1].text(3.2, 0.71 * 5.566 * 0.8, "dynamical-narrowing floor (Entry 4)", fontsize=8, color=INK)
    ax[1].set(xscale="log", yscale="log", xlabel="FWHM, uniform-q control", ylabel="FWHM, profile",
              title="b  No line narrower than a uniform crystal can make")
    ax[1].legend(fontsize=8, frameon=False, loc="upper left")
    for a in ax:
        a.grid(color=GRID, lw=0.6)
        a.set_axisbelow(True)
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
