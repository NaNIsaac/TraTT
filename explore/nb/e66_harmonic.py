"""Entry 66: is the high-frequency excess in Chalise et al. Fig. 2 a second harmonic of the packet?

Full quasi-static TT at the real Si 400 / 17.48 keV parameters (Entry 65), semi-infinite crystal.

Trick: for a packet entering at the surface at t = 0 with v = 1, phi(zeta, t) = S g(t - zeta) depends on
depth only through s = t - zeta (distance behind the front). The Bragg Riccati equation has no other explicit
depth dependence, and below the front the crystal is perfect, so R = r0 (Darwin-Prins root) there. Integrating
dR/ds = +c e^{-i phi} + (i eta_c + c e^{i phi} R) R   (= -dR/dzeta)
from s = 0 (the front) upward gives, at s = t, the surface reflectivity at time t: one ODE per eta gives the
whole R(t). No thickness, no time loop.

Units: length Lambda (extinction length, so alpha = 1), time Lambda / v_L = 1.368 ns, eta = omega.
Weight W(eta) = |Fourier component of R(t) at omega = |eta|| over the window, as in Entry 58 and the paper.
Nonlinearity: K(eta) = (W(eta; S) / S) / (W(eta; S_lin) / S_lin), S_lin = S / 100.
Usage: python e66_harmonic.py [check] [run]
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e65_si400 import params, root  # noqa: E402

P = params(17.48)
T_UNIT = P["Lambda_um"] * 1e-6 / 8433.0          # s per time unit
G_LAMBDA = 2 * np.pi / 1.3577e-10 * P["Lambda_um"] * 1e-6   # G Lambda, Si 400


def to_q(f_ghz):
    return 2 * np.pi * f_ghz * 1e9 * T_UNIT


def g_packet(q0, ncyc):
    """Entry 58 packet: phi = S g, g = exp(-s/tau_d) sin(q0 s) for s > 0. Strain ~ g' starts with a step."""
    tau_d = ncyc * 2 * np.pi / q0
    return lambda s: np.where(s > 0, np.exp(-np.clip(s, 0, None) / tau_d) * np.sin(q0 * s), 0.0)


def surface_R(alpha, kappa, eta, phi, T, per_unit):
    """R(t) at the surface for t in [0, T], semi-infinite crystal; returns t (n,), R (n, n_eta)."""
    eta = np.asarray(eta, float)
    c = 1j * np.pi * alpha
    ie = 1j * (eta + 1j * kappa)
    R = np.array([root(alpha, e, kappa) for e in eta])
    n = int(np.ceil(T * per_unit))
    h = T / n
    s = np.arange(n + 1) * h

    def f(si, r):
        e = np.exp(1j * phi(si))
        return c * np.conj(e) + (ie + c * e * r) * r
    out = np.empty((n + 1, eta.size), complex)
    out[0] = R
    for i in range(n):
        k1 = f(s[i], R)
        k2 = f(s[i] + h / 2, R + h / 2 * k1)
        k3 = f(s[i] + h / 2, R + h / 2 * k2)
        k4 = f(s[i] + h, R + h * k3)
        R = R + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = R
    return s, out


def weight(t, R, eta):
    I = np.abs(R) ** 2
    I = I - I.mean(axis=0)
    w = np.full(t.size, t[1] - t[0])
    w[[0, -1]] /= 2
    return np.abs(np.sum(w[:, None] * I * np.exp(1j * np.abs(eta)[None, :] * t[:, None]), axis=0))


def strain_spec(g, eta, s_max):
    s = np.linspace(0, s_max, 200001)
    eps = np.gradient(g(s), s)
    return np.abs([np.trapezoid(eps * np.exp(1j * o * s), s) for o in np.abs(eta)])


def section_check():
    """Against tratt.core.solve_bragg (slab, back face hidden by kappa = 6) at the Entry 59 case."""
    from tratt.core import required_nz, solve_bragg
    alpha, kappa, S = 8.0, 6.0, 0.3
    q0 = 2 * np.pi * alpha / 0.2
    g = g_packet(q0, 3.0)
    eta = np.array([0.5, 1.0, 2.0]) * q0
    tt = np.linspace(0, 0.6, 61)
    nz = required_nz(alpha, q0 / (2 * np.pi), per_cycle=48, S=S, eta_max=eta.max())
    ref = solve_bragg(alpha, eta, phi=lambda z: S * g(tt[:, None] - z), n_z=nz, kappa=kappa)
    s, R = surface_R(alpha, kappa, eta, lambda x: S * g(x), 0.6, 20000)
    mine = R[np.searchsorted(s, tt - 1e-12)]
    print(f"check vs solve_bragg: max |dR| = {np.max(np.abs(mine - ref)):.2e}, max |R| = {np.max(np.abs(ref)):.2e}")


def section_run():
    out = {}
    f = np.linspace(10, 260, 126)
    eta = to_q(f)
    kappa = 2 * np.pi * P["g"]
    for f0 in (87.0, 113.0):
        q0 = to_q(f0)
        for ncyc in (0.9, 2.0):
            g = g_packet(q0, ncyc)
            S = 1e-3 * G_LAMBDA / q0          # peak strain 1e-3, since max |g'| = q0 at the front
            E = strain_spec(g, eta, 40 * ncyc * 2 * np.pi / q0)
            for T_ps in (140.0, 1400.0):
                T = T_ps * 1e-12 / T_UNIT
                res = {}
                for lab, SS in (("lin", S / 100), ("S", S)):
                    s, R = surface_R(1.0, kappa, eta, lambda x, SS=SS: SS * g(x), T, 64 * (eta.max() + SS * q0) / (2 * np.pi))
                    res[lab] = weight(s, R, eta) / SS
                K = res["S"] / res["lin"]
                key = f"f0={f0:g} ncyc={ncyc:g} T={T_ps:g}ps"
                i0, i2 = np.argmin(abs(f - f0)), np.argmin(abs(f - 2 * f0))
                print(f"{key}: S={S:.2f}  K(f0)={K[i0]:.3f}  K(2f0)={K[i2]:.2f}  K max={K.max():.2f} at {f[np.argmax(K)]:.0f} GHz"
                      f"  |eps(2f0)|/|eps(f0)|={E[i2] / E[i0]:.3f}", flush=True)
                out[key] = dict(S=S, K=K.tolist(), W_S=res["S"].tolist(), W_lin=res["lin"].tolist(), eps=E.tolist())
    json.dump(dict(f=f.tolist(), runs=out), open(HERE / "e66_harmonic.json", "w"))


def section_fig():
    """Like their Fig. 2b/d: strain spectrum vs the reflectivity spectrum with the linear transfer divided out
    (K |eps|), both peak-normalised, so what is left is the strain non-linearity alone."""
    import matplotlib.pyplot as plt
    d = json.load(open(HERE / "e66_harmonic.json"))
    f = np.array(d["f"])
    fig, axs = plt.subplots(2, 2, figsize=(9, 6), sharex=True, constrained_layout=True)
    for i, f0 in enumerate((87, 113)):
        for j, ncyc in enumerate((0.9, 2)):
            ax = axs[i, j]
            r = d["runs"][f"f0={f0} ncyc={ncyc:g} T=140ps"]
            E, K = np.array(r["eps"]), np.array(r["K"])
            ax.plot(f, E / E.max(), "k-", label="strain |ε(f)|")
            ax.plot(f, K * E / (K * E).max(), "--", color="C3", label=f"reflectivity, S = {r['S']:.2f}, linear transfer divided out")
            ax.axvline(2 * f0, color="0.6", lw=0.8, ls=":")
            ax2 = ax.twinx()
            ax2.plot(f, K, color="C0", lw=1)
            ax2.axhline(1, color="C0", lw=0.5, ls=":")
            ax2.set_ylim(0.4, 2.2)
            ax2.tick_params(axis="y", colors="C0")
            ax.set_title(f"{f0} GHz packet, decay {ncyc:g} cycles", fontsize=10)
            if j == 1:
                ax2.set_ylabel("K = non-linear / linear weight", color="C0")
            if j == 0:
                ax.set_ylabel("normalised amplitude")
            if i == 1:
                ax.set_xlabel("frequency (GHz)")
    axs[0, 0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Si 400, 17.48 keV, peak strain 1e-3, 140 ps window, full TT (semi-infinite)", fontsize=11)
    fig.savefig(HERE.parents[1] / "figures" / "si400_harmonic.png", dpi=150)


if __name__ == "__main__":
    for sec in (sys.argv[1:] or ["check", "run"]):
        globals()[f"section_{sec}"]()
