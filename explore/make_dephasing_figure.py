"""Dephasing lifts satellite extinction: chirp (Landau-Zener) and pi phase wall.
Run: python explore/make_dephasing_figure.py

a  n = 1 line at alpha = 2 for spatial chirps beta = 0, 0.1, 0.4: full TT (solid) vs the
   first-Born Fourier model scaled to the same integral as the uniform TT (dotted shape only).
b  plateau height of the chirped satellite vs the Landau-Zener parameter
   Gamma = (pi alpha J_1)^2 / (q0 beta); line: 1 - exp(-2 pi Gamma); kinematic: 2 pi Gamma.
c  integrated satellite with a pi phase wall vs wall depth, relative to no wall (alpha = 0.5, 2);
   kinematic ratio = 1 (Parseval).   (notebook Entries 38-39)
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

from qmod import PhaseDrive, phase_wall, spatial_chirp, uniform  # noqa: E402
from tratt.core import required_nz, rocking  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
HUES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
MARK = ["o", "s", "^", "D", "v"]
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})
AAC, S = 40.0, 1.5
Q0 = 2 * np.pi * AAC
J1 = abs(jv(1, S))
T = 1.5037


def sat(alpha, d, eta, qmax=Q0):
    nz = required_nz(alpha, qmax / (2 * np.pi), S=S, eta_max=eta.max())
    return rocking(alpha, eta, phi=d.phi(T), n_z=nz) - rocking(alpha, eta, n_z=nz)


def main(out="figures/dephasing_extinction.png"):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
    eta = np.arange(0.55 * Q0, 1.45 * Q0, 0.1)
    for beta, c in zip((0.0, 0.1, 0.4), HUES):
        d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=spatial_chirp(Q0, beta) if beta else uniform(Q0))
        I = sat(2.0, d, eta, Q0 * (1 + beta / 2))
        ax[0].plot(eta - Q0, I, color=c, lw=1.2, label=rf"$\beta$={beta:g}  (integral {np.trapezoid(I, eta):.1f})")
    ax[0].set(xlabel=r"$\eta-q_0$", ylabel=r"satellite $I(\eta)$", xlim=(-70, 70),
              title=r"a  Chirp turns the split line into a plateau ($\alpha$=2)")
    ax[0].legend(fontsize=8, frameon=False, loc="upper right")

    for alpha, c, m in zip((0.5, 1.0, 2.0, 4.0), HUES, MARK):
        G, H = [], []
        for beta in (0.05, 0.1, 0.2, 0.4):
            d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=spatial_chirp(Q0, beta))
            I = sat(alpha, d, eta, Q0 * (1 + beta / 2))
            mid = np.abs(eta - Q0) < 0.3 * Q0 * beta / 2
            G.append((np.pi * alpha * J1) ** 2 / (Q0 * beta))
            H.append(np.median(I[mid]))
        ax[1].plot(G, H, m, ms=7, color=c, ls="none", label=rf"$\alpha$={alpha:g}")
    g = np.geomspace(5e-3, 5, 200)
    ax[1].plot(g, 1 - np.exp(-2 * np.pi * g), color=INK, lw=1.3, label=r"Landau–Zener $1-e^{-2\pi\Gamma}$")
    ax[1].plot(g, 2 * np.pi * g, color=INK2, lw=1.0, ls=":", label=r"kinematic $2\pi\Gamma$")
    ax[1].set(xscale="log", ylim=(0, 1.15), xlabel=r"$\Gamma=(\pi\alpha J_1)^2/(q_0\beta)$",
              ylabel="plateau height (full TT)", title="b  Plateau = Landau–Zener transfer, no fit")
    ax[1].legend(fontsize=8, frameon=False, loc="upper left")

    zb = np.linspace(0.0, 1.0, 21)
    for alpha, c, m in zip((0.5, 2.0), HUES, MARK):
        d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=phase_wall(Q0, np.pi, 1.0))
        d0 = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=phase_wall(Q0, 0.0, 1.0))
        t = 1.0 + zb + 0.0037
        nz = required_nz(alpha, AAC, S=S, eta_max=eta.max())
        I0 = rocking(alpha, eta, n_z=nz)[None, :]
        p = np.trapezoid(rocking_cube(alpha, d, eta, t, n_z=nz) - I0, eta, axis=1)
        pu = np.trapezoid(rocking_cube(alpha, d0, eta, t, n_z=nz) - I0, eta, axis=1)
        ax[2].plot(zb, p / pu, m + "-", ms=5, color=c, label=rf"$\alpha$={alpha:g} ($\alpha_1$={alpha*J1:.2f})")
    ax[2].axhline(1.0, color=INK2, lw=1.0, ls=":", label="kinematic (Parseval)")
    ax[2].set(xlabel=r"phase-wall depth $\zeta_b$", ylabel="integrated satellite / no wall",
              title=r"c  A $\pi$ phase wall raises the integral (n=1)")
    ax[2].legend(fontsize=8, frameon=False, loc="upper center")
    for a in ax:
        a.grid(color=GRID, lw=0.6)
        a.set_axisbelow(True)
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
