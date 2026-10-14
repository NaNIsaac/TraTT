"""Entry 65: the Entry 59/64 transfer function for a real crystal, Si 400 symmetric Bragg (Chalise et al. case).

Mapping to the dimensionless Bragg problem of tratt.core.solve_bragg (sigma polarisation, X0h constants of
Entry 20):
    Lambda   = lambda sin(theta_B) / |chi_rh|      (code's Lambda for gamma = sin theta_B)
    g        = kappa / (2 pi alpha) = chi_i0 / |chi_rh|           (Darwin absorption parameter)
    x        = eta / (2 pi alpha) = f / f_D,   f_D = v_L / Lambda  (phonon frequency in Darwin half-widths)
With mu = 2 pi chi_i0 / lambda and kappa = mu L / sin(theta_B), both factors of sin(theta_B) cancel in g.

The crystal factor is the single-mode Fourier weight at the satellite, eta = W, divided by its kinematic
limit |R0| pi alpha S / (2 kappa): H(x, g) = 1 in a thin-absorber kinematic crystal, < 1 where extinction
competes with absorption. It does not depend on the strain spectrum (Entry 59, one-curve collapse).
The heterodyne factor |r0| ~ 1/x is separate and kinematic; H is what remains after it.

Usage: python e65_si400.py [tt] [fig]     tt: full-TT check of the closed form at (g, x) of Si 400
"""
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.core import solve_bragg  # noqa: E402

V_L = 8433.0          # longitudinal sound speed along [100] Si, m/s (Chalise et al. use the same order)
# X0h, Si 400, sigma (Entry 20): E keV -> (theta_B deg, |chi_rh|, chi_i0)
SI400 = {8.0: (34.802, 0.78073e-5, 0.35783e-6), 17.48: (15.142, 0.16049e-5, 0.16059e-7)}
HC = 12.398419843     # keV * Angstrom


def params(E):
    th, xr, xi0 = SI400[E]
    lam = HC / E * 1e-10
    Lam = lam * np.sin(np.radians(th)) / xr
    return dict(E=E, Lambda_um=Lam * 1e6, f_D_GHz=V_L / Lam / 1e9, g=xi0 / xr)


def root(alpha, eta, kappa):
    c = 1j * np.pi * alpha
    r = np.roots([c, 1j * (eta + 1j * kappa), c])
    return r[np.argmin(np.abs(r))]


def H(x, g, alpha=8.0, S=0.02):
    """Closed-form crystal factor at the satellite (Entry 64 single-mode Riccati form)."""
    kappa = g * 2 * np.pi * alpha
    c = 1j * np.pi * alpha
    out = []
    for xi in np.atleast_1d(x):
        eta = W = xi * 2 * np.pi * alpha
        R0 = root(alpha, eta, kappa)
        lam = 1j * (eta + 1j * kappa) + 2 * c * R0
        D = {s: 1j * c * S * (1 - R0 ** 2) / (2 * (lam - 1j * s * W)) for s in (+1, -1)}
        F = np.conj(R0) * D[-1] + R0 * np.conj(D[+1])
        out.append(abs(F) / (abs(R0) * np.pi * alpha * S / (2 * kappa)))
    return np.array(out)


def tt_check(g, xs, S=0.02, nt=32):
    """Full TT, single mode, at eta = W: TT weight / closed form. alpha chosen so kappa = 6 (back face e^-6)."""
    alpha = 6.0 / (2 * np.pi * g)
    kappa = 6.0
    rows = []
    for x in xs:
        W = x * 2 * np.pi * alpha
        t = 2 * np.pi * np.arange(nt) / (W * nt)
        nz = int(24 * (W / (2 * np.pi) + alpha)) + 256

        def phi(zeta):
            return S * np.cos(W * (t[:, None] - zeta))
        R = np.abs(solve_bragg(alpha, np.array([W]), phi=phi, n_z=nz, kappa=kappa)) ** 2
        F_tt = abs(np.mean((R - R.mean(axis=0)) * np.exp(1j * W * t[:, None]), axis=0)[0])
        R0 = root(alpha, W, kappa)
        F_cf = H(x, g, alpha=alpha, S=S)[0] * abs(R0) * np.pi * alpha * S / (2 * kappa)
        rows.append(dict(g=g, x=x, alpha=alpha, n_z=nz, ratio=F_tt / F_cf))
        print(f"  g={g:.4f} x={x:6.1f} (alpha={alpha:.1f}, n_z={nz}): TT / closed form = {F_tt / F_cf:.4f}", flush=True)
    return rows


if __name__ == "__main__":
    P = [params(E) for E in SI400]
    for p in P:
        print(f"Si 400 at {p['E']:g} keV: Lambda = {p['Lambda_um']:.2f} um, f_D = {p['f_D_GHz']:.3f} GHz, g = {p['g']:.4f}")
    fs = np.array([5, 10, 20, 40, 60, 87, 100, 113, 150, 170, 300])
    table = {}
    for p in P:
        h = H(fs / p["f_D_GHz"], p["g"])
        table[p["E"]] = h.tolist()
        print(f"  {p['E']:g} keV  H(f):  " + "  ".join(f"{f}GHz:{v:.3f}" for f, v in zip(fs, h)))
        i20, i100 = list(fs).index(20), list(fs).index(100)
        print(f"     H(20 GHz)/H(100 GHz) = {h[i20] / h[i100]:.3f}")
    out = dict(params=P, f_GHz=fs.tolist(), H=table)
    if "tt" in sys.argv:
        out["tt"] = []
        for p in P:
            out["tt"] += tt_check(p["g"], [3.0, 10.0, 30.0, 100.0, 150.0])
        json.dump(out, open(HERE / "e65_si400.json", "w"), indent=1)
    if "fig" in sys.argv:
        import matplotlib.pyplot as plt
        f = np.geomspace(1, 400, 300)
        fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
        ax.axvspan(20, 170, color="#e4e3df", lw=0, label="Chalise et al. Fig. 2 band")
        for p, c in zip(P, ("#2a78d6", "#eb6834")):
            ax.plot(f, H(f / p["f_D_GHz"], p["g"]), color=c, lw=2.2,
                    label=rf"Si 400, {p['E']:g} keV ($\kappa/2\pi\alpha$ = {p['g']:.3f})")
        ax.axvline(100, color="#52514e", lw=0.8, ls=":")
        ax.set(xscale="log", xlabel="phonon frequency (GHz)", ylabel=r"crystal factor $H$ (weight / kinematic)",
               ylim=(0, 1.05), title="Fourier weight of R(t) relative to the kinematic law")
        ax.legend(frameon=False, fontsize=9, loc="lower right")
        ax.grid(color="#e4e3df", lw=0.6)
        fig.savefig(HERE.parents[1] / "figures" / "si400_transfer.png", dpi=200)
        print("wrote figures/si400_transfer.png")
