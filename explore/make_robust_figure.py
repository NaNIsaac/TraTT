"""Robustness of the transient-only satellite: one control parameter. Run: python explore/make_robust_figure.py

Every non-ideal drive (attenuation, second harmonic, ringing) is reduced to
r = |Theta_1(final)| / max_t |Theta_1(t)| from the accumulated-coupling law, and the
full-TT contrast C = I_max / I_final is plotted against it. Small-angle guide: C = 1/r^2.
Panel b: the essence test - with attenuation, retuning S so that the *depth integral*
of J_1 vanishes restores the contrast although J_1 != 0 over most of the crystal.
Laue, n = 1, alpha_ac = 40; FAST = (alpha 1, T 0.5), MID = (alpha 0.5, T 2) (notebook Entry 35).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from scipy.special import jn_zeros, jv

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "explore" / "nb"))

from e35_robust import AAC, CASES, run  # noqa: E402
from nonideal import GeneralDrive, ramp, ringing_ramp, theta_law  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
mpl.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": "white",
    "savefig.facecolor": "white", "savefig.dpi": 200,
})
Z = jn_zeros(1, 1)[0]


def r_of(d, o):
    th = np.abs(theta_law(d, 1, o["t"]))
    return th[81:].mean() / th[:81].max()


def main(out="figures/robust_contrast.png"):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    fam = {"attenuation": (BLUE, "o"), "second harmonic": (ORANGE, "s"), "ringing (window)": (AQUA, "^")}
    for case, fill in (("FAST", True), ("MID", False)):
        alpha, T = CASES[case]
        drives = {
            "attenuation": [GeneralDrive(S=Z, alpha_ac=AAC, amp=ramp(T), gamma=-np.log(1 - x))
                            for x in (0.01, 0.03, 0.05, 0.1, 0.2, 0.4)],
            "second harmonic": [GeneralDrive(S=Z, alpha_ac=AAC, amp=ramp(T), eps2=x)
                                for x in (0.005, 0.01, 0.02, 0.05, 0.1, 0.2)],
            "ringing (window)": [GeneralDrive(S=Z, alpha_ac=AAC, amp=ringing_ramp(T, R, tau, 3.0))
                                 for R, tau in ((0.05, 0.2), (0.1, 0.2), (0.2, 0.2), (0.1, 1.0), (0.1, 3.0))],
        }
        for name, ds in drives.items():
            c, m = fam[name]
            rs, Cs = [], []
            for d in ds:
                o = run(d, alpha, T, law=False)
                rs.append(r_of(d, o))
                Cs.append(o["C"])
            ax[0].plot(rs, Cs, m, ms=6, color=c, mfc=c if fill else "white", mew=1.2, ls="none",
                       label=f"{name}, {case}")
    r = np.geomspace(1e-3, 1.5, 200)
    ax[0].plot(r, 1 / r**2, color=INK, lw=1.2, label=r"$1/r^2$")
    for C, lab in ((30, "30:1"), (10, "10:1")):
        ax[0].axhline(C, color=INK2, lw=0.8, ls=":")
        ax[0].text(1.1e-3, C * 1.15, lab, fontsize=8, color=INK2)
    ax[0].set(xscale="log", yscale="log", xlim=(1e-3, 1.5), ylim=(0.5, 1e7),
              xlabel=r"$r=|\Theta_1^{\rm final}|\,/\,\max_t|\Theta_1(t)|$ (law)",
              ylabel=r"contrast $I_1^{\max}/I_1^{\rm final}$ (full TT)",
              title="a  All non-idealities collapse onto one parameter")
    ax[0].legend(fontsize=7, frameon=False, loc="upper right", ncol=1)

    # b: essence test
    z = np.linspace(0, 1, 4001)
    losses = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4])
    alpha, T = CASES["FAST"]
    C_bz, C_rt = [], []
    for loss in losses:
        g = -np.log(1 - loss) if loss else 0.0
        Sst = brentq(lambda S: np.trapezoid(jv(1, S * np.exp(-g * z)), z), 3.5, 6.0) if loss else Z
        C_bz.append(run(GeneralDrive(S=Z, alpha_ac=AAC, amp=ramp(T), gamma=g), alpha, T, law=False)["C"])
        C_rt.append(run(GeneralDrive(S=Sst, alpha_ac=AAC, amp=ramp(T), gamma=g), alpha, T, law=False)["C"])
    ax[1].plot(100 * losses, C_bz, "o-", color=BLUE, ms=6, label=r"$S=j_{1,1}$ (local Bessel zero at entrance)")
    ax[1].plot(100 * losses, C_rt, "s-", color=ORANGE, ms=6,
               label=r"$S^*$: $\int_0^1 J_1(S^*e^{-\gamma\zeta})\,d\zeta=0$")
    ax[1].axhline(30, color=INK2, lw=0.8, ls=":")
    ax[1].set(yscale="log", xlabel="acoustic amplitude loss across the crystal (%)",
              ylabel=r"contrast $I_1^{\max}/I_1^{\rm final}$ (full TT)",
              title="b  What must vanish is the depth integral, not $J_1$")
    ax[1].legend(fontsize=8, frameon=False, loc="center right")
    for a in ax:
        a.grid(color=GRID, lw=0.6, which="major")
        a.set_axisbelow(True)
    fig.savefig(ROOT / out)
    print("wrote", out)


if __name__ == "__main__":
    main()
