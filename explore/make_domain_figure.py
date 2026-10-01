"""Two acoustic domains: coherent cancellation at the centre, incoherent sum in the integral.
Run: python explore/make_domain_figure.py

Transducer: S2 = 5.331 (J_1 min) for one transit, then S1 = 1.841 (J_1 max). The domain wall
moves in at v_s; at wall depth r* = 0.373 the two contributions to Theta_1 cancel.
a  centre I_1(t): both domains vs each domain alone (Laue, alpha = 2).
b  line shape at the cancellation vs the S1 domain alone: the centre is dark, the
   satellite is not (window integral 17.6 vs 16.6) - Parseval, notebook Entry 37.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import jv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "explore" / "nb"))

from memory_drive import HistoryDrive  # noqa: E402
from tratt.core import required_nz, rocking  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})
AAC, ALPHA = 40.0, 2.0
S1 = minimize_scalar(lambda s: -jv(1, s), bounds=(1, 3), method="bounded").x
S2 = minimize_scalar(lambda s: jv(1, s), bounds=(4, 7), method="bounded").x
R = jv(1, S2) / (jv(1, S2) - jv(1, S1))


def drive(a2, a1):
    return HistoryDrive(S=1.0, alpha_ac=AAC, env=lambda s: np.where(s < 1.0, a2, a1))


def main(out="figures/domain_cancellation.png"):
    q = 2 * np.pi * AAC
    nz = required_nz(ALPHA, AAC, S=S2, eta_max=1.5 * q)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    t = np.linspace(0.9, 2.05, 231) + 0.0037
    for (a2, a1), c, lab in (((S2, S1), INK, "both domains"), ((0.0, S1), BLUE, r"$S_1$ domain alone"),
                             ((S2, 0.0), ORANGE, r"$S_2$ domain alone")):
        I = rocking_cube(ALPHA, drive(a2, a1), np.array([q]), t, n_z=nz)[:, 0]
        ax[0].plot(t - 1, I, color=c, lw=1.8 if c == INK else 1.3, label=lab)
    ax[0].axvline(R, color=INK2, lw=0.8, ls="--")
    ax[0].text(R + 0.01, 0.5, rf"wall at $r^*$={R:.3f}", fontsize=8, color=INK2)
    ax[0].set(xlabel=r"wall depth $\zeta_b=\tilde t-1$", ylabel=r"$I(q,\tilde t)$ (centre)",
              xlim=(-0.1, 1.05), title=r"a  Each domain diffracts, together they cancel ($n=1$, centre)")
    ax[0].legend(fontsize=8, frameon=True, framealpha=0.9, edgecolor="none", loc="lower right")

    eta = np.linspace(q - 40, q + 40, 1601)
    ts = 1.0 + R + 0.0037
    Ib = rocking(ALPHA, eta, phi=drive(S2, S1).phi(ts), n_z=nz)
    Ia = rocking(ALPHA, eta, phi=drive(0.0, S1).phi(ts), n_z=nz)
    ax[1].plot(eta - q, Ib, color=INK, lw=1.6, label=f"both (window integral {np.trapezoid(Ib, eta):.1f})")
    ax[1].plot(eta - q, Ia, color=BLUE, lw=1.2, label=rf"$S_1$ alone ({np.trapezoid(Ia, eta):.1f})")
    ax[1].set(xlabel=r"$\eta-q$", ylabel=r"$I(\eta)$ at the cancellation",
              title="b  The centre is dark, the satellite is not")
    ax[1].legend(fontsize=8, frameon=False, loc="upper right")
    for a in ax:
        a.grid(color=GRID, lw=0.6)
        a.set_axisbelow(True)
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
