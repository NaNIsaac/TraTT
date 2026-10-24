"""Entry 77 (T2c): anomalous absorption chi_ih in Bragg, and whether the crystal factor H needs it.

solve_bragg now takes kappa_h (chi_h'' in kappa units, kappa_h / kappa = chi_ih / chi_i0); for a centrosymmetric
crystal both couplings become c = i pi alpha + kappa_h / 2. The first-order closed form of Entries 59/65/75 goes
through with complex c: R0 from c R0^2 + i (eta + i kappa) R0 + c = 0, Lambda = i (eta + i kappa) + 2 c R0,
D_s = i c S (1 - R0^2) / (2 (Lambda - i s W)), F = conj(R0) D_- + R0 conj(D_+), H = |F| / (|R0| pi alpha S / 2 kappa).
The detuning eta = +-W now matters (Borrmann asymmetry), so H is given on both sides.

Pre-registered rule (notebook Entry 77): Si 400, X0h ratio kappa_h / kappa = 0.93; if max over 20-150 GHz of
|H(kappa_h) / H(0) - 1| exceeds 1 % at g = 0.0100 (17.48 keV, either sign of eta), the closed form needs chi_ih as a
second parameter. Check: the extended closed form against full TT (single mode, as e65 tt) to 1 %.
Usage: python e77_chi_ih.py [run] [tt] [fig]
"""

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e65_si400 import params  # noqa: E402
from e75_window_cf import packet_specs, setup  # noqa: E402
from tratt.core import solve_bragg  # noqa: E402

OUT = HERE / "e77_chi_ih.json"
RATIO = 0.93  # chi_ih / chi_i0, Si 400 (X0h: 0.925 at 8 keV, 0.931 at 17.48 keV; Entry 20 script)
F = np.array([5, 10, 20, 30, 40, 60, 87, 100, 150, 200, 300], float)
BAND = (F >= 20) & (F <= 150)


def root(alpha, eta, kappa, kappa_h=0.0):
    c = 1j * np.pi * alpha + 0.5 * kappa_h
    r = np.roots([c, 1j * (eta + 1j * kappa), c])
    return r[np.argmin(np.abs(r))], c


def H_steady(x, g, ratio, sign, alpha=8.0, S=0.02):
    """Crystal factor at detuning eta = sign * x * 2 pi alpha, mode frequency W = |eta|."""
    kappa = g * 2 * np.pi * alpha
    out = []
    for xi in np.atleast_1d(x):
        W = xi * 2 * np.pi * alpha
        eta = sign * W
        R0, c = root(alpha, eta, kappa, ratio * kappa)
        lam = 1j * (eta + 1j * kappa) + 2 * c * R0
        D = {s: 1j * c * S * (1 - R0**2) / (2 * (lam - 1j * s * W)) for s in (+1, -1)}
        Fw = np.conj(R0) * D[-1] + R0 * np.conj(D[+1])
        out.append(abs(Fw) / (abs(R0) * np.pi * alpha * S / (2 * kappa)))
    return np.array(out)


def H_window(E, f, T_ns, ratio, sign):
    """CF-B of Entry 75 with complex c, Entry 66 packet, window T_ns; same normalisation as e75 (|G_T| |R0| pi S / kappa)."""
    st = setup(E)
    kappa, S = st["kappa"], st["S"]
    T = T_ns * 1e-9 / st["t_unit"]
    out = []
    for fj in np.atleast_1d(f):
        W = 2 * np.pi * fj * 1e9 * st["t_unit"]
        eta = sign * W
        R0, c = root(1.0, eta, kappa, ratio * kappa)
        lam = 1j * (eta + 1j * kappa) + 2 * c * R0
        G, L = packet_specs(st["q0"], st["tau"], eta, lam, T)
        Ft = {Wk: 1j * c * S * (1 - R0**2) * (G(Wk) - np.exp((lam + 1j * Wk) * T) * L) / (lam + 1j * Wk) for Wk in (W, -W)}
        Iw = np.conj(R0) * Ft[W] + R0 * np.conj(Ft[-W])
        out.append(abs(Iw) / (abs(G(W)) * abs(R0) * np.pi * S / kappa))
    return np.array(out)


def section_run():
    out = dict(ratio=RATIO, f_GHz=F.tolist(), cases={})
    for E in (17.48, 8.0):
        p = params(E)
        g = p["g"]
        x = F / p["f_D_GHz"]
        case = dict(g=g, f_D=p["f_D_GHz"])
        print(f"Si 400 {E:g} keV, g = {g:.4f}, kappa_h/kappa = {RATIO}")
        for sign, lab in ((+1, "eta>0"), (-1, "eta<0")):
            H0 = H_steady(x, g, 0.0, sign)
            H1 = H_steady(x, g, RATIO, sign)
            r = H1 / H0
            case[f"steady_{lab}"] = dict(H0=H0.tolist(), H1=H1.tolist(), ratio=r.tolist())
            print(f"  steady {lab}: H(kh)/H(0) = " + " ".join(f"{fj:g}:{v:.4f}" for fj, v in zip(F, r, strict=True))
                  + f"   max|r-1| in band {np.max(np.abs(r[BAND] - 1)) * 100:.2f} %, outside {np.max(np.abs(r[~BAND] - 1)) * 100:.2f} %")
            for T_ns in (14.0, 52.0):
                Hw0 = H_window(E, F, T_ns, 0.0, sign)
                Hw1 = H_window(E, F, T_ns, RATIO, sign)
                rw = Hw1 / Hw0
                case[f"T{T_ns:g}_{lab}"] = dict(H0=Hw0.tolist(), H1=Hw1.tolist(), ratio=rw.tolist())
                print(f"  window {T_ns:g} ns {lab}: H(kh)/H(0) = " + " ".join(f"{fj:g}:{v:.4f}" for fj, v in zip(F, rw, strict=True))
                      + f"   max|r-1| in band {np.max(np.abs(rw[BAND] - 1)) * 100:.2f} %")
        asym = np.array(case["steady_eta<0"]["H1"]) / np.array(case["steady_eta>0"]["H1"])
        case["asym_steady"] = asym.tolist()
        print("  steady H(eta<0)/H(eta>0) with kh: " + " ".join(f"{fj:g}:{v:.4f}" for fj, v in zip(F, asym, strict=True)))
        out["cases"][f"E{E:g}"] = case
    c17 = out["cases"]["E17.48"]
    worst = {k: float(np.max(np.abs(np.array(c17[k]["ratio"])[BAND] - 1)) * 100) for k in c17 if isinstance(c17[k], dict)}
    out["worst_band_17keV_pct"] = worst
    wmax = max(worst.values())
    out["verdict"] = "second parameter needed" if wmax > 1.0 else "one parameter suffices (<= 1 % in band)"
    print(f"verdict at g = {c17['g']:.4f}: max |H(kh)/H(0) - 1| over 20-150 GHz, steady and windows 14/52 ns, both signs: "
          + ", ".join(f"{k} {v:.2f} %" for k, v in worst.items()) + f" -> max {wmax:.2f} % -> {out['verdict']}")
    json.dump(out, open(OUT, "w"), indent=1)


def section_tt(signs=(+1, -1)):
    """Single mode at (g, x), kappa_h = RATIO kappa, both signs: TT weight / extended closed form. Settings
    from Entry 65 sec. 3: S = 0.002 (linear, nu < 0.03) and 48 RK4 steps per cycle (the 24-step, S = 0.02
    run of e65 is not converged at x >= 30 and must not be quoted)."""
    g = params(17.48)["g"]
    alpha = 6.0 / (2 * np.pi * g)
    kappa = 6.0
    S, nt = 0.002, 32
    out = json.load(open(OUT)) if OUT.exists() else {}
    rows = [r for r in out.get("tt", []) if r["sign"] not in signs]
    for sign in signs:
        for x in (3.0, 10.0, 30.0, 60.0):
            t0 = time.time()
            W = x * 2 * np.pi * alpha
            eta = sign * W
            t = 2 * np.pi * np.arange(nt) / (W * nt)
            nz = int(48 * (W / (2 * np.pi) + alpha)) + 256

            def phi(zeta):
                return S * np.cos(W * (t[:, None] - zeta))  # noqa: B023

            R = np.abs(solve_bragg(alpha, np.array([eta]), phi=phi, n_z=nz, kappa=kappa, kappa_h=RATIO * kappa)) ** 2
            F_tt = abs(np.mean((R - R.mean(axis=0)) * np.exp(1j * W * t[:, None]), axis=0)[0])
            R0, c = root(alpha, eta, kappa, RATIO * kappa)
            F_cf = H_steady(x, g, RATIO, sign, alpha=alpha, S=S)[0] * abs(R0) * np.pi * alpha * S / (2 * kappa)
            # same run without kappa_h, for the size of the effect in TT itself
            R_0 = np.abs(solve_bragg(alpha, np.array([eta]), phi=phi, n_z=nz, kappa=kappa)) ** 2
            F_tt0 = abs(np.mean((R_0 - R_0.mean(axis=0)) * np.exp(1j * W * t[:, None]), axis=0)[0])
            F_cf0 = H_steady(x, g, 0.0, sign, alpha=alpha, S=S)[0] * abs(root(alpha, eta, kappa)[0]) * np.pi * alpha * S / (2 * kappa)
            rows.append(dict(sign=sign, x=x, nz=nz, tt_over_cf=F_tt / F_cf, tt_ratio=F_tt / F_tt0, cf_ratio=F_cf / F_cf0))
            print(f"  g={g:.4f} x={x:6.1f} eta{'>' if sign > 0 else '<'}0 (n_z={nz}): TT / closed form = {F_tt / F_cf:.4f};"
                  f"  effect of kappa_h: TT {F_tt / F_tt0:.4f}, closed form {F_cf / F_cf0:.4f} ({time.time() - t0:.0f} s)", flush=True)
    out["tt"] = rows
    out["tt_worst_pct"] = float(max(abs(r["tt_over_cf"] - 1) for r in rows) * 100)
    print(f"TT vs extended closed form: worst {out['tt_worst_pct']:.2f} %")
    json.dump(out, open(OUT, "w"), indent=1)


def section_conv60():
    """Step convergence at x = 60, eta > 0: 48 -> 96 RK4 steps per cycle, with and without kappa_h, so the
    1 % residual of the tt check at x = 60 can be attributed (numerics / S) rather than to the extension."""
    g = params(17.48)["g"]
    alpha = 6.0 / (2 * np.pi * g)
    kappa, S, nt, x = 6.0, 0.002, 32, 60.0
    W = x * 2 * np.pi * alpha
    t = 2 * np.pi * np.arange(nt) / (W * nt)
    out = json.load(open(OUT))
    res = {}
    for per in (48, 96):
        nz = int(per * (W / (2 * np.pi) + alpha)) + 256
        for kh in (0.0, RATIO):
            t0 = time.time()
            R = np.abs(solve_bragg(alpha, np.array([W]), phi=lambda z: S * np.cos(W * (t[:, None] - z)), n_z=nz, kappa=kappa, kappa_h=kh * kappa)) ** 2
            F_tt = abs(np.mean((R - R.mean(axis=0)) * np.exp(1j * W * t[:, None]), axis=0)[0])
            R0, _ = root(alpha, W, kappa, kh * kappa)
            F_cf = H_steady(x, g, kh, +1, alpha=alpha, S=S)[0] * abs(R0) * np.pi * alpha * S / (2 * kappa)
            res[f"per{per}_kh{kh:g}"] = F_tt / F_cf
            print(f"  x=60 eta>0 {per} steps/cycle (n_z={nz}) kappa_h/kappa={kh:g}: TT / closed form = {F_tt / F_cf:.4f} ({time.time() - t0:.0f} s)", flush=True)
    out["conv60"] = res
    json.dump(out, open(OUT, "w"), indent=1)


def section_ttp():
    section_tt(signs=(+1,))


def section_ttm():
    section_tt(signs=(-1,))


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    f = np.array(d["f_GHz"])
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    ax = axs[0]
    c = d["cases"]["E17.48"]
    ax.plot(f, c["steady_eta>0"]["H0"], "k-", lw=1.5, label="$\\chi_{ih} = 0$ (Entry 65)")
    ax.plot(f, c["steady_eta>0"]["H1"], "C0-", lw=1.5, label="$\\chi_{ih} = 0.93\\,\\chi_{i0}$, $\\eta > 0$")
    ax.plot(f, c["steady_eta<0"]["H1"], "C3-", lw=1.5, label="$\\chi_{ih} = 0.93\\,\\chi_{i0}$, $\\eta < 0$ (Borrmann edge)")
    ax.axvspan(20, 150, color="#e4e3df", lw=0, zorder=0)
    ax.set(xscale="log", xlabel="phonon frequency (GHz)", ylabel="crystal factor $H$ (steady)", title="Si 400, 17.48 keV", ylim=(0, 1.05))
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.grid(color="#e4e3df", lw=0.6)
    ax = axs[1]
    for E, ls in (("17.48", "-"), ("8", "--")):
        c = d["cases"][f"E{E}"]
        for lab, col in (("eta>0", "C0"), ("eta<0", "C3")):
            ax.plot(f, 100 * (np.array(c[f"steady_{lab}"]["ratio"]) - 1), ls, color=col, lw=1.5, label=f"{E} keV, {lab}")
    ax.axhspan(-1, 1, color="#e4e3df", lw=0, zorder=0)
    ax.axvspan(20, 150, color="#f2f1ee", lw=0, zorder=0)
    ax.axhline(0, color="0.5", lw=0.6)
    ax.set(xscale="log", xlabel="phonon frequency (GHz)", ylabel="$H(\\chi_{ih}) / H(0) - 1$ (%)", title="change of $H$ with anomalous absorption")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(color="#e4e3df", lw=0.6)
    fig.savefig(HERE.parents[1] / "figures" / "paper" / "e77_chi_ih.png", dpi=200)
    print("wrote figures/paper/e77_chi_ih.png")


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "tt", "fig"]:
        globals()[f"section_{sec}"]()
