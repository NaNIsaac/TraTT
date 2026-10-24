"""Entry 82 (T2k): what beam divergence does to the crystal factor H(f; T_w).

A ray at detuning eta' = eta + delta rings at its own frequency ~eta' with damping kappa after the packet has passed,
so the Fourier weight at the nominal omega = |eta| collects only the rays within the ringing linewidth. Beam average:
    W_conv(omega) = | sum_k p_k I~(eta + delta_k; omega, T_w) |     (intensities add, the oscillating parts add as phasors)
with p a Gaussian in delta of width sigma_eta = 2 pi sigma (sigma in Darwin half-widths). TT: every ray is a full
semi-infinite TT run (Entry 75 machinery, vectorised over eta). Closed form: CF-B of Entry 75 at detuning eta + delta_k
and Fourier frequency omega = |eta|.

Pre-registered rule (notebook Entry 82): |W_TT,conv / W_CF,conv - 1| > 1 % at any (f, sigma, T_w) -> "divergence term
needed"; otherwise "closed form retained". Reported: H_conv(f; T_w, sigma) normalised at 100 GHz against T_w, T_flat /
T_90 where they exist, and H_conv / H_single at 300 ns against the expectation kappa sqrt(pi/2) / sigma_eta.
Usage: python e82_divergence.py [run] [fig]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e65_si400 import SI400, params, root  # noqa: E402
from e75_window_cf import integrate_packet, packet_specs, setup  # noqa: E402

OUT = HERE / "e82_divergence.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e82_divergence.png"
F = np.array([20.0, 40.0, 100.0])
T_NS = (1.4, 14.0, 52.0, 300.0)
SIGMAS = (0.0, 0.3, 1.0, 3.0)  # Darwin half-widths


def delta_grid(sigma_eta, kappa, du=0.025):
    """Non-uniform grid delta = s sinh(u): spacing kappa/80 at the centre, ~0.025 |delta| in the tails, to +-6 sigma.
    (The first run used du = 0.1 and failed the 0.2 % grid-convergence check by 2 %; du = 0.025 passes, see log.)"""
    if sigma_eta == 0.0:
        return np.array([0.0]), np.array([1.0])
    s = kappa / 2
    U = np.arcsinh(6 * sigma_eta / s)
    u = np.arange(-U, U + du / 2, du)
    d = s * np.sinh(u)
    w = np.exp(-0.5 * (d / sigma_eta) ** 2)
    dd = np.gradient(d)
    p = w * dd
    return d, p / p.sum()


def fourier_complex(st, run, T, omega):
    """Complex Fourier component at omega (per eta) of |R|^2 - |R0|^2 over [0, T]; e75.tt_weight without the abs."""
    s, R, r1, r2, Lam = run["s"], run["R"], run["r1"], run["r2"], run["Lam"]
    I = np.abs(R) ** 2 - np.abs(r1) ** 2
    if T <= run["s_p"]:
        n = int(round(T / run["h"]))
        w = np.full(n + 1, run["h"])
        w[[0, -1]] /= 2
        return np.sum(w[:, None] * I[: n + 1] * np.exp(1j * omega[None, :] * s[: n + 1, None]), axis=0), n * run["h"]
    w = np.full(s.size, run["h"])
    w[[0, -1]] /= 2
    tot = np.sum(w[:, None] * I * np.exp(1j * omega[None, :] * s[:, None]), axis=0)
    u_p = (R[-1] - r1) / (R[-1] - r2)
    assert np.all(np.abs(u_p) < 0.5)
    K = int(np.ceil(np.log(1e-15) / np.log(np.abs(u_p).max()))) + 1
    k = np.arange(1, K + 1)
    d = (r1 - r2)[None, :] * u_p[None, :] ** k[:, None]
    sp, dT = run["s_p"], T - run["s_p"]

    def fint(a, om):
        b = a + 1j * om
        return np.exp(1j * om * sp) * (np.exp(b * dT) - 1) / b

    lin = np.zeros(omega.size, complex)
    for j in range(omega.size):
        a = k * Lam[j]
        lin[j] = np.conj(r1[j]) * np.sum(d[:, j] * fint(a, omega[j])) + r1[j] * np.sum(np.conj(d[:, j]) * fint(np.conj(a), omega[j]))
        aa = (k[:, None] * Lam[j] + k[None, :] * np.conj(Lam[j])).ravel()
        dd = (d[:, j][:, None] * np.conj(d[:, j])[None, :]).ravel()
        lin[j] += np.sum(dd * fint(aa, omega[j]))
    return tot + lin, T


def cf_complex(st, eta_p, omega, T):
    """CF-B (Entry 75) at detuning eta_p, Fourier frequency omega, window T; complex, per element."""
    kappa, S = st["kappa"], st["S"]
    c = 1j * np.pi
    out = np.empty(eta_p.size, complex)
    for i, (e, W) in enumerate(zip(eta_p, omega, strict=True)):
        R0 = root(1.0, e, kappa)
        lam = 1j * (e + 1j * kappa) + 2 * c * R0
        G, L = packet_specs(st["q0"], st["tau"], e, lam, T)
        Ft = {Wk: 1j * c * S * (1 - R0**2) * (G(Wk) - np.exp((lam + 1j * Wk) * T) * L) / (lam + 1j * Wk) for Wk in (W, -W)}
        out[i] = np.conj(R0) * Ft[W] + R0 * np.conj(Ft[-W])
    return out


def norm(st, omega, T):
    """Kinematic law |G_T(omega)| |R0(omega)| pi S / kappa at the nominal detuning (so that H_single is the Entry 75 H)."""
    kappa, S = st["kappa"], st["S"]
    out = []
    for W in omega:
        R0 = root(1.0, W, kappa)
        G, _ = packet_specs(st["q0"], st["tau"], W, 0.0, T)
        out.append(abs(G(W)) * abs(R0) * np.pi * S / kappa)
    return np.array(out)


def section_run():
    t0 = time.time()
    base = setup(17.48, s_lin=1000.0)
    kappa = base["kappa"]
    eta_nom = 2 * np.pi * F * 1e9 * base["t_unit"]
    out = dict(f_GHz=F.tolist(), T_ns=T_NS, sigmas=SIGMAS, kappa=kappa, cases={}, worst=0.0)
    for sigma in SIGMAS:
        sig_eta = 2 * np.pi * sigma
        d, p = delta_grid(sig_eta, kappa)
        eta_all = (eta_nom[:, None] + d[None, :]).ravel()
        om_all = np.repeat(eta_nom, d.size)
        st = dict(base, eta=eta_all)
        run = integrate_packet(st)
        rows = []
        for T_ns in T_NS:
            T = T_ns * 1e-9 / base["t_unit"]
            Itt, Ta = fourier_complex(st, run, T, om_all)
            Icf = cf_complex(st, eta_all, om_all, Ta)
            Wtt = np.abs((Itt.reshape(F.size, -1) * p[None, :]).sum(axis=1))
            Wcf = np.abs((Icf.reshape(F.size, -1) * p[None, :]).sum(axis=1))
            hn = norm(base, eta_nom, Ta)
            rows.append(dict(T_ns=T_ns, H_tt=(Wtt / hn).tolist(), H_cf=(Wcf / hn).tolist(), ratio=(Wtt / Wcf).tolist()))
            out["worst"] = max(out["worst"], float(np.max(np.abs(Wtt / Wcf - 1))))
        out["cases"][f"s{sigma:g}"] = dict(sigma=sigma, sigma_eta=sig_eta, n_delta=int(d.size), rows=rows)
        print(f"sigma = {sigma:g} Darwin half-widths (sigma_eta = {sig_eta:.2f}, {d.size} rays, RK4 {run['s'].size - 1} steps):", flush=True)
        for r in rows:
            print(f"   T_w = {r['T_ns']:5g} ns: H_conv TT " + " ".join(f"{f:g}:{v:.4f}" for f, v in zip(F, r["H_tt"], strict=True))
                  + "  | TT/CF " + " ".join(f"{v:.4f}" for v in r["ratio"]) + f"  | H(20)/H(100) = {r['H_tt'][0] / r['H_tt'][2]:.3f}")
        if sigma > 0:
            exp_flat = kappa * np.sqrt(np.pi / 2) / sig_eta
            print(f"   expectation at 300 ns, flat: kappa sqrt(pi/2)/sigma_eta = {exp_flat:.4f};  H_conv(300 ns) = "
                  + " ".join(f"{v:.4f}" for v in rows[-1]["H_tt"]) + f"  ({time.time() - t0:.0f} s)")
            out["cases"][f"s{sigma:g}"]["expect_flat"] = float(exp_flat)
    out["verdict"] = "closed form retained" if out["worst"] < 0.01 else "divergence term needed"
    print(f"rule: worst |W_TT,conv / W_CF,conv - 1| = {out['worst'] * 100:.3f} % -> {out['verdict']}")
    # convergence: du halved at sigma = 1, 14 ns
    sig_eta = 2 * np.pi
    for du in (0.1, 0.05, 0.025, 0.0125):
        d, p = delta_grid(sig_eta, kappa, du=du)
        st = dict(base, eta=eta_nom[1] + d)
        run = integrate_packet(st)
        Itt, Ta = fourier_complex(st, run, 14e-9 / base["t_unit"], np.full(d.size, eta_nom[1]))
        W = abs((Itt * p).sum())
        print(f"convergence sigma = 1, 40 GHz, 14 ns, du = {du}: {d.size} rays, W_conv = {W:.6e}")
        out[f"conv_du{du}"] = float(W)
    # angle and energy scale of the ringing width, both energies
    scale = {}
    for E in (17.48, 8.0):
        th, xr, _ = SI400[E]
        g = params(E)["g"]
        thD = xr / np.sin(np.radians(2 * th))  # Darwin half-width, symmetric Bragg, sigma
        dth = g * thD
        dE = dth / np.tan(np.radians(th)) * E * 1e6  # meV
        scale[f"E{E:g}"] = dict(darwin_half_urad=thD * 1e6, ringing_half_urad=dth * 1e6, ringing_meV=dE, g=g)
        print(f"E = {E:g} keV: Darwin half-width {thD * 1e6:.2f} urad; ringing width kappa = {g:.4f} half-widths = {dth * 1e6:.3f} urad = {dE:.1f} meV")
    out["scale"] = scale
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT.name} ({time.time() - t0:.0f} s)")


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    ax = axs[0]
    for sig, c in zip(d["sigmas"], ("k", "C0", "C1", "C3"), strict=True):
        rows = d["cases"][f"s{sig:g}"]["rows"]
        T = [r["T_ns"] for r in rows]
        ax.plot(T, [r["H_tt"][0] / r["H_tt"][2] for r in rows], "o-", color=c, label=f"$\\sigma$ = {sig:g} Darwin half-widths")
    ax.axhline(0.511, color="0.6", lw=0.8, ls="--")
    ax.set(xscale="log", xlabel="window $T_w$ (ns)", ylabel="$H(20\\,\\mathrm{GHz})/H(100\\,\\mathrm{GHz})$", title="the crystal factor's shape against beam divergence", ylim=(0.45, 1.1))
    ax.legend(fontsize=8, frameon=False)
    ax.grid(color="#e4e3df", lw=0.6)
    ax = axs[1]
    for sig, c in zip(d["sigmas"], ("k", "C0", "C1", "C3"), strict=True):
        rows = d["cases"][f"s{sig:g}"]["rows"]
        T = [r["T_ns"] for r in rows]
        ax.plot(T, [r["H_tt"][0] for r in rows], "o-", color=c, label=f"$\\sigma$ = {sig:g}, 20 GHz")
        ax.plot(T, [r["H_tt"][2] for r in rows], "s--", color=c, mfc="none", label=f"$\\sigma$ = {sig:g}, 100 GHz")
    ax.set(xscale="log", yscale="log", xlabel="window $T_w$ (ns)", ylabel="$H_{conv}$ (weight / kinematic law)", title="absolute $H_{conv}$: divergence sets the coherence time")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    ax.grid(color="#e4e3df", lw=0.6)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
