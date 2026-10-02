"""Entry 59: first-order (linear-in-S) closed form for the Fourier weight of R(t) in absorbing Bragg,
compared with the full TT weights of e58_chalise_ref.py (sections abs, darwin). No TT run here.

Linearising the Riccati equation of tratt.core.solve_bragg around the semi-infinite static root R0
(c R0^2 + i eta_c R0 + c = 0, |R0| < 1, eta_c = eta + i kappa, c = i pi alpha), with phi = S g(t - zeta):
    d delta / d zeta = i c phi (1 - R0^2) - Lambda delta,   Lambda = i eta_c + 2 c R0,
    int delta(0, t) e^{i W t} dt = i c S (1 - R0^2) G(W) / (Lambda + i W),   G(W) = int g(tau) e^{i W tau} dtau.
The weight measured in e58 (mean over the window T of (R - <R>) e^{i W t}, W = |eta|) is then
    F = [ conj(R0) A(W) + R0 conj(A(-W)) ] / T,   A(W) = i c S (1 - R0^2) G(W) / (Lambda + i W).
Kinematic tail: |R0| ~ pi alpha / eta, Lambda - i eta ~ -kappa  ->  W ~ |R0| |G| (pi alpha S / kappa) / T.
Usage: python e59_transfer.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e58_chalise_ref import packet  # noqa: E402

T = 0.6


def root(alpha, eta, kappa):
    c = 1j * np.pi * alpha
    out = []
    for e in np.atleast_1d(eta):
        r = np.roots([c, 1j * (e + 1j * kappa), c])
        out.append(r[np.argmin(np.abs(r))])
    return np.array(out)


def G(g, tau_d, w):
    tau = np.linspace(0, 40 * tau_d, 40001)
    gg = g(tau)
    return np.array([np.trapezoid(gg * np.exp(1j * x * tau), tau) for x in np.atleast_1d(w)])


def predict(alpha, kappa, S, q0, eta):
    c = 1j * np.pi * alpha
    g, tau_d = packet(S, q0)
    R0 = root(alpha, eta, kappa)
    W = np.abs(eta)
    lam = 1j * (eta + 1j * kappa) + 2 * c * R0
    A = lambda w: 1j * c * S * (1 - R0 ** 2) * G(g, tau_d, w) / (lam + 1j * w)  # noqa: E731
    return np.abs(np.conj(R0) * A(W) + R0 * np.conj(A(-W))) / T, R0


def compare(lab, alpha, kappa, S, q0, r):
    eta = np.array(r["eta"])
    Wn = np.array(r["W"])
    Wp, R0 = predict(alpha, kappa, S, q0, eta)
    rel = Wn / Wp
    wd = 2 * np.pi * alpha
    sel = [np.argmin(abs(eta / wd - x)) for x in (0.4, 0.7, 1.0, 1.5, 2, 3, 5, 10, 20, 40, 60)]
    sel = sorted(set(i for i in sel if eta[i] / wd >= 0.3))
    print(f" {lab}: W_TT / W_pred at eta/(2 pi alpha) = " +
          "  ".join(f"{eta[i] / wd:.1f}:{rel[i]:.3f}" for i in sel))
    big = Wn > 0.02 * Wn.max()
    print(f"   where W > 2% of max: median {np.median(rel[big]):.4f}, range [{rel[big].min():.3f}, "
          f"{rel[big].max():.3f}];  kinematic-only law (|R0| |G| pi alpha S / kappa / T) off by up to "
          f"{np.max(np.abs(np.abs(R0[big]) * np.abs(G(*packet(S, q0), eta[big])) * np.pi * alpha * S / kappa / T / Wn[big] - 1)):.2f}")
    return rel


if __name__ == "__main__":
    alpha, kappa = 8.0, 6.0
    A = json.load(open(HERE / "e58_abs.json"))
    for S in (0.02, 0.05):
        compare(f"R_D=0.05 S={S}", alpha, kappa, S, 2 * np.pi * alpha / 0.05, A[f"S{S}"])
    f = HERE / "e58_darwin.json"
    if f.exists():
        D = json.load(open(f))
        for RD in (0.2, 0.5):
            compare(f"R_D={RD} S=0.02", alpha, kappa, 0.02, D[f"RD{RD}"]["q0"], D[f"RD{RD}"])


def figure(out=HERE.parents[1] / "figures" / "chalise_transfer.png"):
    """W from full TT (points), the first-order closed form (lines) and the kinematic tail law (dashed)."""
    import matplotlib.pyplot as plt
    alpha, kappa, S = 8.0, 6.0, 0.02
    wd = 2 * np.pi * alpha
    D = json.load(open(HERE / "e58_darwin.json"))
    A = json.load(open(HERE / "e58_abs.json"))
    cases = [(0.05, A["S0.02"], 2 * np.pi * alpha / 0.05), (0.2, D["RD0.2"], D["RD0.2"]["q0"]),
             (0.5, D["RD0.5"], D["RD0.5"]["q0"])]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    for (RD, r, q0), c in zip(cases, ("#2a78d6", "#eb6834", "#1baf7a")):
        eta = np.array(r["eta"])
        Wn = np.array(r["W"])
        Wp, R0 = predict(alpha, kappa, S, q0, eta)
        Wk = np.abs(R0) * np.abs(G(*packet(S, q0), eta)) * np.pi * alpha * S / kappa / T
        x = eta / wd
        ax[0].plot(x, Wn, "o", ms=3.5, color=c, mec="white", mew=0.4, label=rf"TT, $R_D$={RD}")
        ax[0].plot(x, Wp, "-", color=c, lw=1.4)
        ax[0].plot(x, Wk, "--", color=c, lw=1.0)
        ax[1].plot(x, Wn / Wk, "o", ms=3.5, color=c, mec="white", mew=0.4, label=rf"$R_D$={RD}")
        ax[1].plot(x, Wp / Wk, "-", color=c, lw=1.4)
    ax[0].set(xscale="log", yscale="log", xlabel=r"$\eta/(2\pi\alpha)$ = phonon frequency / Darwin half-width",
              ylabel=r"Fourier weight $W$ of $R(t)$",
              title="a  TT (points), first-order closed form (lines), kinematic law (dashed)")
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].axhline(1, color="#52514e", lw=0.8, ls=":")
    ax[1].set(xscale="log", yscale="log", xlabel=r"$\eta/(2\pi\alpha)$", ylabel="W / kinematic law",
              title=r"b  W / kinematic law: one curve ($\kappa/2\pi\alpha$=0.12)")
    ax[1].legend(frameon=False, fontsize=8)
    for a in ax:
        a.grid(color="#e4e3df", lw=0.6)
        a.spines[["top", "right"]].set_visible(False)
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__" and "fig" in sys.argv:
    figure()


def weak_absorption_check():
    """TT check of the closed form at kappa / (2 pi alpha) = 0.024 (alpha = 40, kappa = 6, R_D = 0.2)."""
    from e58_chalise_ref import run
    alpha, kappa, S = 40.0, 6.0, 0.02
    wd = 2 * np.pi * alpha
    q0 = wd / 0.2
    eta = np.linspace(0.3 * wd, 3 * q0, 60)
    r = run(alpha, S, q0, eta, kappa=kappa)
    Wp, R0 = predict(alpha, kappa, S, q0, eta)
    rel = r["W"] / Wp
    big = r["W"] > 0.02 * r["W"].max()
    print(f" weak absorption (kappa/2pi alpha = {kappa / wd:.3f}), n_z = {r['nz']}: W_TT / W_pred median "
          f"{np.median(rel[big]):.4f}, range [{rel[big].min():.3f}, {rel[big].max():.3f}]")
    json.dump({k: np.asarray(v).tolist() for k, v in r.items()} | {"pred": Wp.tolist()},
              open(HERE / "e59_weak.json", "w"))


if __name__ == "__main__" and "weak" in sys.argv:
    weak_absorption_check()
