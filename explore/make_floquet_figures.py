"""Figures for notebook Entries 44-50. Run: python explore/make_floquet_figures.py

floquet_structure.png      a reduced-zone bands, b gaps vs R, c Bloch-mixing points, d peak trajectories
reduced_model_convergence.png  E_pop vs retained forward sidebands, and k_req vs R
time_regimes.png           I(eta, t) for resolved / strong overlap / dark main line + window traces
anomaly_summary.png        E_ind and E_G3 vs R; darkest strong and brightest near-zero orders
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
from scipy.special import jv  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "explore" / "nb"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(NB))
from e46_bands import bands  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
HUES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
mpl.rcParams.update({"font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": "white",
                     "savefig.facecolor": "white", "savefig.dpi": 200})


def grid(ax):
    for a in np.atleast_1d(ax).flat:
        a.grid(color=GRID, lw=0.6)
        a.set_axisbelow(True)


def floquet_structure():
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.4), constrained_layout=True)
    aac, S = 5.0, 1.5
    q = 2 * np.pi * aac
    eta = np.linspace(-0.5 * q, 2.5 * q, 1201)
    for R, c in zip((0.05, 0.5, 1.0), HUES):
        lam, gap, mix, _ = bands(R * aac, aac, S, eta)
        ax[0].plot(eta / q, lam[:, 0] / q, ",", color=c)
        ax[0].plot(eta / q, lam[:, 1] / q, ",", color=c)
        ax[0].plot([], [], "o", color=c, label=f"R={R:g}")
    ax[0].set(xlabel=r"$\eta/q$", ylabel=r"quasi-propagation constant $\lambda/q$ (reduced zone)",
              title="a  Two exact Floquet branches (S=1.5)")
    ax[0].legend(fontsize=8, frameon=False, markerscale=1)
    Rs = np.linspace(0.02, 1.5, 30)
    for n, c in zip(range(4), HUES):
        g = []
        for R in Rs:
            _, gp, _, _ = bands(R * aac, aac, S, np.array([n * q]))
            g.append(gp[0])
        ax[1].plot(Rs, np.array(g) / q, "o", ms=3.5, color=c, label=f"n={n}")
        ax[1].plot(Rs, np.minimum(2 * np.pi * Rs * aac * abs(jv(n, S)), np.inf) / q, color=c, lw=0.8)
    ax[1].axhline(0.5, color=INK2, lw=0.8, ls=":")
    ax[1].set(xlabel="R", ylabel="gap at $\\eta=nq$ / q", ylim=(0, 0.6),
              title=r"b  Gaps vs independent $2\pi\alpha|J_n|$ (lines)")
    ax[1].legend(fontsize=8, frameon=False)
    for S2, c in zip((0.5, 1.0, 2.3), HUES):
        d = np.load(NB / f"e48_traj_S{S2}.npz", allow_pickle=True)
        for R, cr in zip(d["Rs"], d["cross"]):
            for x in cr:
                if x < 0.49:  # crossings at q/2 are sorting artefacts (notebook Entry 48)
                    ax[2].plot(R, x, "o", ms=3.5, color=c)
        ax[2].plot([], [], "o", color=c, label=f"S={S2:g}")
        for R, pk in zip(d["Rs"], d["peaks"]):
            for x in pk:
                if 0 <= x <= 2.6:
                    ax[3].plot(R, x, "o", ms=2.5, color=c)
        ax[3].plot([], [], "o", color=c, label=f"S={S2:g}")
    ax[2].set(xlabel=r"R ($\alpha_{ac}$=2)", ylabel=r"Bloch-mixing point $\eta^*/q$",
              title="c  Bloch-mixing points move and annihilate")
    ax[2].legend(fontsize=8, frameon=False)
    for k in range(3):
        ax[3].axhline(k, color=INK2, lw=0.6, ls=":")
    ax[3].set(xlabel=r"R ($\alpha_{ac}$=2)", ylabel=r"spectral maxima $\eta/q$ (phase-averaged)",
              title="d  Peaks leave the channel positions")
    ax[3].legend(fontsize=8, frameon=False)
    grid(ax)
    fig.savefig(ROOT / "figures" / "floquet_structure.png")


def convergence():
    rows = json.load(open(NB / "e44b_sidebands.json"))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for r in rows:
        if r["S"] == 1.5 and r["aac"] == 5.0:
            ax[0].plot(range(len(r["errs"])), np.maximum(r["errs"], 1e-4), "o-", ms=3,
                       color=plt.cm.viridis(min(r["R"] / 3, 1)), label=f"R={r['R']:g}")
    ax[0].axhline(0.02, color=INK2, lw=0.8, ls=":")
    ax[0].set(yscale="log", xlabel="forward-beam sidebands kept, |k| <= kmax",
              ylabel=r"$E_{pop}$ vs converged", title=r"a  Convergence in Floquet modes (S=1.5, $\alpha_{ac}$=5)")
    ax[0].legend(fontsize=7, frameon=False, ncol=2)
    for (S, aac), c in zip(((1.5, 5.0), (1.5, 2.0), (4.0, 5.0), (4.0, 2.0)), HUES):
        pts = [(r["R"], r["kreq"] if r["kreq"] is not None else 9) for r in rows if r["S"] == S and r["aac"] == aac]
        ax[1].plot(*zip(*pts), "o-", color=c, ms=4, label=rf"S={S:g}, $\alpha_{{ac}}$={aac:g}")
    Rf = np.linspace(0, 3, 50)
    ax[1].plot(Rf, 4 * Rf, color=INK2, lw=0.8, ls="--", label=r"$4R$")
    ax[1].set(xlabel="R", ylabel=r"$k_{req}$ ($E_{pop}<0.02$; 9 = not reached)", ylim=(-0.3, 9.5),
              title="b  Sidebands needed grow with R, saturate at the Bessel width")
    ax[1].legend(fontsize=8, frameon=False)
    grid(ax)
    fig.savefig(ROOT / "figures" / "reduced_model_convergence.png")


def time_regimes():
    d = np.load(NB / "e49_time.npy", allow_pickle=True).item()
    names = ["1 resolved (R 0.1)", "3 strong overlap (R 1)", "5 dark main line (R 1.05, S 0.8)"]
    fig, ax = plt.subplots(2, 3, figsize=(16, 7.5), constrained_layout=True)
    from tratt.core import required_nz
    from tratt.cube import rocking_cube
    from tratt.drive import Drive
    from e49_time import CASES
    aac = 5.0
    q = 2 * np.pi * aac
    for j, name in enumerate(names):
        alpha, S = CASES[name]
        eta = np.arange(-2.5 * q, 2.5 * q, 0.3)
        t = np.linspace(0.0, 1.2, 97) + 0.0037
        I = rocking_cube(alpha, Drive(S=S, alpha_ac=aac), eta, t, n_z=required_nz(alpha, aac, S=S, eta_max=2.5 * q))
        ax[0, j].imshow(I, origin="lower", aspect="auto", cmap="Blues", extent=[eta[0] / q, eta[-1] / q, t[0], t[-1]],
                        vmin=0, vmax=np.percentile(I, 99.5))
        ax[0, j].set(title=name, xlabel=r"$\eta/q$", ylabel=r"$\tilde t$")
        tt, _, res = d[name]
        for (n, p, pind, c), col in zip(res, HUES):
            ax[1, j].plot(tt, p, color=col, lw=1.2, label=f"n={n} full")
            ax[1, j].plot(tt, pind, color=col, lw=0.8, ls="--")
        ax[1, j].set(xlabel=r"$\tilde t$", ylabel="window integral (dashed: independent)")
        ax[1, j].legend(fontsize=7, frameon=False)
    grid(ax[1])
    fig.savefig(ROOT / "figures" / "time_regimes.png")


def anomaly():
    rows = json.load(open(NB / "e50_anomaly.json"))
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.4), constrained_layout=True)
    R = np.array([r["R"] for r in rows])
    ax[0].plot(R, [r["E_ind"] for r in rows], "o", ms=3, color=HUES[0], label="independent channels")
    ax[0].plot(R, [r["E_g3"] for r in rows], "s", ms=3, mfc="white", color=HUES[1], label="third-order couplings only")
    ax[0].set(xscale="log", yscale="log", xlabel="R", ylabel=r"$E_{pop}$ (phase-averaged)",
              title="a  The naive third-order model is not a better global model")
    ax[0].legend(fontsize=8, frameon=False)
    dk = [r for r in rows if np.isfinite(r["dark"])]
    ax[1].plot([r["R"] for r in dk], [r["dark"] for r in dk], "o", ms=3, color=HUES[2],
               label=r"darkest order with $|J_n|>0.2$ ($p_n/p_n^{G3}$)")
    br = [r for r in rows if np.isfinite(r["bright"])]
    ax[1].plot([r["R"] for r in br], [r["bright"] for r in br], "^", ms=3, color=HUES[3],
               label=r"brightest order with $|J_n|<0.05$ ($p_n/p_{max}$)")
    ax[1].set(xscale="log", yscale="log", xlabel="R", title="b  Dark / bright anomalies appear only at R >~ 1")
    ax[1].legend(fontsize=8, frameon=False, loc="lower left")
    grid(ax)
    fig.savefig(ROOT / "figures" / "anomaly_summary.png")


if __name__ == "__main__":
    for f in (convergence, anomaly, floquet_structure, time_regimes):
        f()
        print("done", f.__name__, flush=True)
