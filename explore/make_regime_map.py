"""Regime map of the filled-depth law. Run: python explore/make_regime_map.py

Axes are the two independent dimensionless knobs:
  x = Lambda_ac / Lambda_P = alpha / alpha_ac   -> how clean the mapping is (resolution)
  y = alpha_n = L |J_n(S)| / Lambda_P           -> how many satellite Pendellosung turns fit in the transit
Colour: max over t in [0, 1] of |I(nq, t) - sin^2(pi alpha_n t)|, Laue, n = 1, S = 1.5.
In the stroboscopic model v_s enters only through t_tr: at fixed material the drive
frequency moves x alone, and the clock tau_n = Lambda_P / (v_s |J_n|) does not depend on it.
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
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})


def main(out="figures/regime_map.png"):
    S, n = 1.5, 1
    J = abs(jv(n, S))
    xs = np.geomspace(0.0125, 0.4, 7)
    ys = np.geomspace(0.25, 4.0, 7)
    t = np.linspace(0.0, 1.0, 121) + 0.0037
    E = np.empty((ys.size, xs.size))
    for i, y in enumerate(ys):
        alpha = y / J
        for j, x in enumerate(xs):
            aac = alpha / x
            d = Drive(S=S, alpha_ac=aac)
            nz = int(max(1500, 40 * max(aac, alpha)))
            I = rocking_cube(alpha, d, np.array([n * d.q]), t, n_z=nz)[:, 0]
            E[i, j] = np.max(np.abs(I - np.sin(np.pi * y * t) ** 2))
            print(f"x={x:.4f} alpha_1={y:.3f} alpha={alpha:.3f} aac={aac:.1f} err={E[i, j]:.3f}")

    fig, ax = plt.subplots(figsize=(6.6, 5.0), constrained_layout=True)
    im = ax.pcolormesh(xs, ys, E, cmap="Blues", shading="nearest", vmin=0, vmax=0.5)
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            ax.text(x, y, f"{E[i, j]:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if E[i, j] > 0.3 else INK)
    ax.set(xscale="log", yscale="log",
           xlabel=r"$\Lambda_{\rm ac}/\Lambda_P=\alpha/\alpha_{\rm ac}$  (set by drive frequency)",
           ylabel=r"$\alpha_1=L|J_1(S)|/\Lambda_P$  (turns in transit)",
           title=r"Error of $\sin^2(\pi\alpha_1\tilde t)$, Laue $n=1$, $S=1.5$")
    fig.colorbar(im, ax=ax, label=r"max$_{\tilde t}\,|I-\sin^2|$")
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
