"""Time-resolved rocking curve near satellite n = 1: full TT vs the closed three-mode model.
Run: python explore/make_map_figure.py

Laue, S = 1.5, alpha = 3, alpha_ac = 60 (alpha / alpha_ac = 0.05, notebook Entry 19).
Third panel: centre cut and the sin^2 law. Colour is a single-hue sequential map
(magnitude), shared scale on both maps.
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
sys.path.insert(0, str(ROOT / "explore" / "nb"))

from e14_threemode import model  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

INK, INK2 = "#0b0b0b", "#52514e"
BLUE, ORANGE = "#2a78d6", "#eb6834"
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})


def main(out="figures/poster_map.png"):
    S, alpha, n, aac = 1.5, 3.0, 1, 60.0
    d = Drive(S=S, alpha_ac=aac)
    q = d.q
    an = alpha * abs(jv(n, S))
    deta = np.linspace(-40, 40, 321)
    eta = n * q + deta
    t = np.linspace(0.02, 1.0, 99) + 0.0037
    I = rocking_cube(alpha, d, eta, t, n_z=int(60 * aac))
    M = np.array([model(alpha, S, q, n, eta, x) for x in t])
    print(f"max |model - TT| / max TT = {np.max(np.abs(M - I)) / I.max():.3f}")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True,
                           gridspec_kw=dict(width_ratios=[1, 1, 0.8]))
    kw = dict(origin="lower", aspect="auto", cmap="Blues", vmin=0, vmax=I.max(),
              extent=[deta[0], deta[-1], t[0], t[-1]])
    im = ax[0].imshow(I, **kw)
    ax[1].imshow(M, **kw)
    ax[0].set(title="a  Full Takagi–Taupin", xlabel=r"$\eta-q$", ylabel=r"$\tilde t=t/t_{\rm tr}$")
    ax[1].set(title="b  Three-mode two-layer model (closed form)", xlabel=r"$\eta-q$")
    fig.colorbar(im, ax=ax[:2], shrink=0.85, label=r"$I(\eta,\tilde t)$")
    k0 = np.argmin(np.abs(deta))
    tf = np.linspace(0, 1, 400)
    ax[2].plot(tf, np.sin(np.pi * an * tf) ** 2, color=INK, lw=1.4, label=r"$\sin^2(\pi\alpha_1\tilde t)$")
    ax[2].plot(t, I[:, k0], "o", ms=3.5, color=BLUE, mec="white", mew=0.5, label="full TT")
    ax[2].plot(t, M[:, k0], "s", ms=3, color=ORANGE, mec="white", mew=0.5, label="three-mode")
    ax[2].set(title=r"c  Centre cut $\eta=q$", xlabel=r"$\tilde t$", ylabel=r"$I(q,\tilde t)$",
              ylim=(-0.03, 1.08))
    for s in ("top", "right"):
        ax[2].spines[s].set_visible(False)
    ax[2].legend(fontsize=8, frameon=False, loc="lower right")
    fig.suptitle(rf"Laue, satellite $n=1$, $S$={S}, $\alpha$={alpha:g}, $\alpha_{{\rm ac}}$={aac:g}"
                 rf" ($\alpha_1$={an:.2f}): the satellite narrows and blinks as the front fills the crystal",
                 fontsize=10)
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
