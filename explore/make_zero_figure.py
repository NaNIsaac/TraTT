"""Bessel-zero resurrection and the renormalised zero. Run: python explore/make_zero_figure.py

a  channel population at S = j_{n,1} (exact Floquet) vs alpha_ac for n = 0, 1, 2 and alpha = 0.5, 1, 2;
   lines: third-order theory sin^2|g_eff|, g_eff = pi alpha (R/2)^2 C_n(j).
b  collapse: sqrt(P) / (pi alpha) vs (R/2)^2 |C_n|  (identity line).
c  renormalised zero: S* where the channel is darkest vs R; line: j + (R/2)^2 C_n / J'_{-n}(j).
(notebook Entry 45)
"""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.special import jvp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "explore" / "nb"))
from e45_zero import ZEROS, C_n, channel_population  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
HUES = {0: "#2a78d6", 1: "#eb6834", 2: "#1baf7a"}
MARK = {0.5: "o", 1.0: "s", 2.0: "^"}
mpl.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": "white",
                     "savefig.facecolor": "white", "savefig.dpi": 200})
# S* values measured in Entry 45b (n, R, S*_numeric)
SSTAR = [(0, 0.20, 2.392924), (0, 0.333, 2.370842), (0, 0.40, 2.356748), (0, 0.667, 2.268143), (0, 0.50, 2.327234),
         (0, 1.00, 2.066517), (1, 0.20, 3.822626), (1, 0.333, 3.805810), (1, 0.40, 3.795193), (1, 0.667, 3.727901)]


def main(out="figures/bessel_zero_resurrection.png"):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
    aacs = np.array([40, 20, 10, 5, 3, 2, 1.5, 1.0])
    for n, j in ZEROS.items():
        Cn = C_n(n, j)
        for alpha in (0.5, 1.0, 2.0):
            P = np.array([channel_population(alpha, a, j, n)[0] for a in aacs])
            R = alpha / aacs
            g = np.pi * alpha * (R / 2) ** 2 * Cn
            ax[0].plot(aacs, P, MARK[alpha], ms=6, color=HUES[n], ls="none",
                       label=rf"$J_{n}=0$, $\alpha$={alpha:g}" if True else None)
            af = np.geomspace(1, 40, 100)
            ax[0].plot(af, np.sin(np.pi * alpha * (alpha / af / 2) ** 2 * abs(Cn)) ** 2, color=HUES[n], lw=0.8)
            ax[1].plot((R / 2) ** 2 * abs(Cn), np.arcsin(np.sqrt(P)) / (np.pi * alpha), MARK[alpha], ms=6,
                       color=HUES[n], ls="none")
    x = np.geomspace(1e-5, 0.5, 50)
    ax[1].plot(x, x, color=INK, lw=1.0, label="third order")
    ax[0].set(xscale="log", yscale="log", xlabel=r"$\alpha_{ac}$", ylabel=r"channel population at $J_n(S)=0$",
              title="a  A dark order comes back as channels overlap")
    ax[0].legend(fontsize=7, frameon=False, ncol=3, loc="lower left")
    ax[1].set(xscale="log", yscale="log", xlabel=r"$(R/2)^2|C_n(j)|$", ylabel=r"$\arcsin\sqrt{P}\,/\,\pi\alpha$",
              title=r"b  Collapse: $g_{\rm eff}=\pi\alpha\,(R/2)^2 C_n$ (no fit)")
    ax[1].legend(fontsize=8, frameon=False, loc="upper left")
    Rf = np.linspace(0, 1.05, 100)
    for n in (0, 1):
        j = ZEROS[n]
        dJ = jvp(-n, j) if n else jvp(0, j)
        ax[2].plot(Rf, j + (Rf / 2) ** 2 * C_n(n, j) / dJ - j, color=HUES[n], lw=1.2, label=rf"theory, $n={n}$")
        pts = [(R, s) for nn, R, s in SSTAR if nn == n]
        ax[2].plot([p[0] for p in pts], [p[1] - j for p in pts], "o", color=HUES[n], ms=6, label=rf"full, $n={n}$")
    ax[2].set(xlabel=r"$R=\alpha/\alpha_{ac}$", ylabel=r"$S^*-j_{n,1}$",
              title="c  The darkest amplitude moves: renormalised zero")
    ax[2].legend(fontsize=8, frameon=False, loc="lower left")
    for a in ax:
        a.grid(color=GRID, lw=0.6)
        a.set_axisbelow(True)
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
