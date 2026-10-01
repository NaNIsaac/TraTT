"""Dark loci of orders n = 0, 1, 2 in the (R, S) plane. Run: python explore/make_darkmap_figure.py

D_n = p_n / p_n^independent (phase-averaged window integral over the independent-channel value),
Laue, alpha_ac = 2 (top) and 5 (bottom); notebook Entry 47.  Overlays: Bessel zeros of J_n
(dashed), the third-order renormalised zero J_{-n}(S) = (R/2)^2 C_n(S) (dotted), and
R |J_n(S)| = 1 (solid thin).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jn_zeros, jv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "explore" / "nb"))
from e45_zero import C_n  # noqa: E402

INK, INK2 = "#0b0b0b", "#52514e"
mpl.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": "white",
                     "savefig.facecolor": "white", "savefig.dpi": 200})


def main(out="figures/darkmap.png"):
    fig, ax = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True, sharex=True, sharey=True)
    Sf = np.linspace(0.2, 6.0, 300)
    for r, aac in enumerate((2.0, 5.0)):
        d = np.load(ROOT / "explore" / "nb" / f"e47b_darkmap_q{aac}.npz")
        Rs, Ss, D = d["Rs"], d["Ss"], d["D"]
        for j, n in enumerate((0, 1, 2)):
            a = ax[r, j]
            im = a.pcolormesh(Ss, Rs, np.log10(np.clip(D[:, :, j], 1e-3, 10)), cmap="RdBu", vmin=-2, vmax=1,
                              shading="nearest")
            for z in jn_zeros(n, 2):
                a.axvline(z, color=INK, lw=0.8, ls="--")
            # third-order renormalised zero near the first Bessel zero
            z0 = jn_zeros(n, 1)[0]
            Sz = np.linspace(z0 - 1.2, z0 + 0.3, 400)
            g = [(jv(-n, s) / C_n(n, s)) for s in Sz]
            Rz = np.sqrt(np.clip(4 * np.array(g), 0, None))
            ok = (Rz > 0) & (Rz < 1.5)
            a.plot(Sz[ok], Rz[ok], color=INK, lw=1.0, ls=":")
            Rj = 1 / np.abs(jv(n, Sf))
            a.plot(Sf[Rj < 1.5], Rj[Rj < 1.5], color=INK2, lw=0.7)
            a.set(title=rf"$n={n}$, $\alpha_{{ac}}$={aac:g}", xlim=(0.2, 6), ylim=(0.05, 1.5))
            if r == 1:
                a.set_xlabel("S")
            if j == 0:
                a.set_ylabel(r"$R=\alpha/\alpha_{ac}=\Lambda_{ac}/\Lambda_P$")
    fig.colorbar(im, ax=ax, shrink=0.7, label=r"$\log_{10} D_n$  (red: darker, blue: brighter than independent)")
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
