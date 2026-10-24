"""Entry 70: does the detuning of maximum signal coincide with the detuning of maximum Fisher information?

Directive question B2, one case only: Si 400 symmetric Bragg at 17.48 keV (g = 0.0100) and, as the
crystal-factor control, 8 keV (g = 0.0458); the Entry 66 packet (87 GHz, decay 0.9 cycles, peak strain
1e-3) in the semi-infinite retarded-time solver of Entry 66. The measurement is R(t) at a fixed detuning
eta (one angle), sampled uniformly over a window T, with additive Gaussian noise of constant sigma.

Signal: the Fourier weight W(eta) at omega = |eta|, as in Chalise et al. and Entries 58-66.
Information: for the strain parameters theta = (S, f0, tau_d), the per-detuning Fisher information
    I_ij(eta) = (1/sigma^2) sum_t dR/dtheta_i dR/dtheta_j,   R = |R(t)|^2,
with central finite differences (relative step 1 %, checked at 0.5 %). Reported per parameter as the
diagonal I_ii(eta) and as the marginal precision 1/[I^-1]_ii(eta) (the other two parameters unknown).

Pre-registered rule (written before the run). eta_W = argmax W(eta); eta_I(theta_i) = argmax I_ii(eta).
  - For S the trivial expectation is eta_I = eta_W (R is linear in S, so dR/dS is the signal itself).
    "Different" is claimed only if |f_I(S) - f_W| exceeds both one Darwin half-width (0.73 GHz) and 10 %
    of f0 (8.7 GHz); otherwise the answer to B2 for the amplitude is "they coincide".
  - For f0 and tau_d a flank maximum is expected from Fourier theory alone and is not a finding. What is
    reported is whether the crystal moves it: the shift between 17.48 keV and 8 keV, and between the full
    model and a kinematic reference (same packet, H = 1: the reflectivity replaced by the heterodyne
    |r0| * displacement spectrum |eps^|/omega, Entry 59 sec. 1). A crystal effect is claimed only if the 17.48 keV flank maximum differs
    from the 8 keV one by more than 10 % of f0.
Window: 140 ps (their Fig. 2a) primary, 1400 ps check.
Usage: python e70_fisher.py [run] [fig]
"""

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e65_si400 import params, root  # noqa: E402
from e66_harmonic import surface_R  # noqa: E402

V_L = 8433.0
F_MIN = float(os.environ.get("E70_FMIN", "10"))  # post-hoc rerun uses 30 (Entry 70)
F_GRID = np.arange(F_MIN, 261.0, 1.0)  # GHz, detuning scan
TAG = "" if F_MIN == 10 else f"_f{F_MIN:g}"
OUT = HERE / f"e70_fisher{TAG}.json"
FIG = HERE.parents[1] / "figures" / "paper" / f"e70_fisher{TAG}.png"
THETA0 = dict(
    S=None, f0=87.0, tau_d=None
)  # S and tau_d filled per energy (strain 1e-3, 0.9 cycles)


def units(E):
    p = params(E)
    t_unit = p["Lambda_um"] * 1e-6 / V_L  # s per time unit
    g_lambda = 2 * np.pi / 1.3577e-10 * p["Lambda_um"] * 1e-6  # G Lambda, Si 400
    return p, t_unit, g_lambda


def to_q(f_ghz, t_unit):
    return 2 * np.pi * f_ghz * 1e9 * t_unit


def packet(q0, tau_d):
    return lambda s: np.where(s > 0, np.exp(-np.clip(s, 0, None) / tau_d) * np.sin(q0 * s), 0.0)


def intensity(E, theta, T_ps, per_cycle=64, per_unit=None):
    """|R(t)|^2 on the common time grid for parameters theta = (S, f0 GHz, tau_d in time units)."""
    p, t_unit, _ = units(E)
    S, f0, tau_d = theta
    q0 = to_q(f0, t_unit)
    eta = to_q(F_GRID, t_unit)
    kappa = 2 * np.pi * p["g"]
    T = T_ps * 1e-12 / t_unit
    per_unit = per_unit or 1.05 * per_cycle * (eta.max() + S * q0) / (2 * np.pi)
    g = packet(q0, tau_d)
    t, R = surface_R(1.0, kappa, eta, lambda s: S * g(s), T, per_unit)
    return t, np.abs(R) ** 2, eta


def fourier_weight(t, inten, eta):
    inten = inten - inten.mean(axis=0)
    w = np.full(t.size, t[1] - t[0])
    w[[0, -1]] /= 2
    return np.abs(
        np.sum(w[:, None] * inten * np.exp(1j * np.abs(eta)[None, :] * t[:, None]), axis=0)
    )


def fisher(E, theta, T_ps, rel=0.01):
    """Per-detuning Fisher matrix (sigma = 1) by central differences; returns t, I0, W, Iij (3,3,n_eta)."""
    t, I0, eta = intensity(E, theta, T_ps)
    per_unit = (t.size - 1) / t[-1]  # same time grid for every finite-difference run
    W = fourier_weight(t, I0, eta)
    D = []
    for i in range(3):
        th_p, th_m = list(theta), list(theta)
        h = rel * theta[i]
        th_p[i] += h
        th_m[i] -= h
        _, Ip, _ = intensity(E, th_p, T_ps, per_unit=per_unit)
        _, Im, _ = intensity(E, th_m, T_ps, per_unit=per_unit)
        D.append((Ip - Im) / (2 * h))
    D = np.array(D)  # (3, nt, n_eta)
    Iij = np.einsum("itk,jtk->ijk", D, D)
    return t, W, Iij, eta


def kinematic_reference(E, theta, T_ps):
    """H = 1 reference: oscillation = 2 |r0| * (first-order kinematic satellite), i.e. the weight is
    |r0(eta)| times the strain spectrum of the packet at omega = eta, and its parameter derivatives."""
    p, t_unit, _ = units(E)
    eta = to_q(F_GRID, t_unit)
    kappa = 2 * np.pi * p["g"]
    r0 = np.abs([root(1.0, e, kappa) for e in eta])

    def spec(th):
        S, f0, tau_d = th
        q0 = to_q(f0, t_unit)
        s = np.linspace(0, T_ps * 1e-12 / t_unit, 40001)
        phi = S * packet(q0, tau_d)(s)
        eps = np.gradient(phi, s)
        return (
            np.abs([np.trapezoid(eps * np.exp(1j * o * s), s) for o in eta]) / eta
        )  # |phi^| = |eps^| / omega

    W = r0 * spec(theta)
    D = []
    for i in range(3):
        th_p, th_m = list(theta), list(theta)
        h = 0.01 * theta[i]
        th_p[i] += h
        th_m[i] -= h
        D.append(r0 * (spec(th_p) - spec(th_m)) / (2 * h))
    D = np.array(D)
    return W, D**2  # diagonal information proxy (signal amplitude derivative squared)


def section_run():
    out = {}
    for E in (17.48, 8.0):
        p, t_unit, g_lambda = units(E)
        q0 = to_q(87.0, t_unit)
        theta = (1e-3 * g_lambda / q0, 87.0, 0.9 * 2 * np.pi / q0)
        for T_ps in (140.0, 1400.0):
            _, W, Iij, eta = fisher(E, theta, T_ps)
            Iinv = np.array(
                [np.linalg.pinv(Iij[:, :, k]) for k in range(eta.size)]
            )  # (n_eta, 3, 3)
            prec = 1.0 / np.array([Iinv[k].diagonal() for k in range(eta.size)]).T  # (3, n_eta)
            names = ("S", "f0", "tau_d")
            rec = dict(
                E=E,
                g=p["g"],
                theta=list(theta),
                T_ps=T_ps,
                f=F_GRID.tolist(),
                W=W.tolist(),
                I_diag=[Iij[i, i].tolist() for i in range(3)],
                prec=prec.tolist(),
            )
            fW = F_GRID[np.argmax(W)]
            print(f"E={E} keV g={p['g']:.4f} S={theta[0]:.3f} T={T_ps:g} ps: f_W = {fW:.0f} GHz")
            for i, nme in enumerate(names):
                fI = F_GRID[np.argmax(Iij[i, i])]
                fP = F_GRID[np.argmax(prec[i])]
                print(
                    f"   {nme:6s} argmax I_ii = {fI:.0f} GHz, argmax marginal precision = {fP:.0f} GHz"
                )
                rec[f"fI_{nme}"] = float(fI)
                rec[f"fP_{nme}"] = float(fP)
            rec["fW"] = float(fW)
            if T_ps == 140.0:
                Wk, Dk = kinematic_reference(E, theta, T_ps)
                rec["kin_W"] = Wk.tolist()
                rec["kin_D2"] = Dk.tolist()
                print(
                    f"   kinematic H=1 reference: f_W = {F_GRID[np.argmax(Wk)]:.0f} GHz; argmax D^2: "
                    + ", ".join(f"{n} {F_GRID[np.argmax(Dk[i])]:.0f}" for i, n in enumerate(names))
                )
            out[f"E{E}_T{T_ps:g}"] = rec
        # step check at 140 ps
        _, _, I2, _ = fisher(E, theta, 140.0, rel=0.005)
        I1 = np.array(out[f"E{E}_T140"]["I_diag"])
        d = np.max(np.abs(np.array([I2[i, i] for i in range(3)]) / I1 - 1))
        print(f"   finite-difference step 1 % vs 0.5 %: max relative change of I_ii = {d:.2e}")
        out[f"E{E}_T140"]["step_check"] = float(d)
    json.dump(out, open(OUT, "w"))


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for ax, E in zip(axs, (17.48, 8.0), strict=True):
        r = d[f"E{E}_T140"]
        f = np.array(r["f"])
        W = np.array(r["W"])
        ax.plot(f, W / W.max(), "k-", lw=2, label="Fourier weight W (signal)")
        for i, (nme, c) in enumerate(zip(("S", "f0", "tau_d"), ("C0", "C3", "C2"), strict=True)):
            info = np.array(r["I_diag"][i])
            P = np.array(r["prec"][i])
            ax.plot(f, info / info.max(), "-", color=c, lw=1.2, label=f"Fisher I_ii, {nme}")
            ax.plot(f, P / P.max(), ":", color=c, lw=1.2, label=f"marginal precision, {nme}")
        ax.axvline(87, color="0.6", lw=0.8, ls="--")
        ax.set(
            xlabel="detuning as phonon frequency (GHz)",
            ylabel="normalised",
            title=f"Si 400, {E:g} keV, g = {r['g']:.4f}, 140 ps window",
        )
        ax.grid(color="#e4e3df", lw=0.6)
    axs[0].legend(fontsize=7, loc="upper right")
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


def section_diag():
    """Is the low-frequency maximum of W the second-order (S^2) slow term plus window leakage?
    Compare W at S with W at S/100 scaled by 100 (first order only), 140 ps window, 17.48 keV."""
    p, t_unit, g_lambda = units(17.48)
    q0 = to_q(87.0, t_unit)
    S = 1e-3 * g_lambda / q0
    tau_d = 0.9 * 2 * np.pi / q0
    for lab, SS in (("S", S), ("S/100 x 100", S / 100)):
        t, I0, eta = intensity(17.48, (SS, 87.0, tau_d), 140.0)
        W = fourier_weight(t, I0, eta) * (S / SS)
        i87 = np.argmin(abs(F_GRID - 87))
        print(
            f"  {lab:12s}: f_W = {F_GRID[np.argmax(W)]:.0f} GHz, W(f_min)/W(87) = {W[0] / W[i87]:.2f}, "
            f"cycles in window at f_min = {F_GRID[0] * 1e9 * 140e-12:.1f}"
        )


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
