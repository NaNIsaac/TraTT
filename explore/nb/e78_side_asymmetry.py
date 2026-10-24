"""Entry 78 (T2f): the +-eta side asymmetry of the crystal factor with chi_ih, tested with a packet in full TT.

Entry 77 (closed form, single mode): with chi_ih the two detunings +-eta of one phonon frequency get different
crystal factors, and A = H(-eta) / H(+eta) cancels the strain spectrum. Here: full semi-infinite TT with the
Entry 66 packet (Entry 75 machinery, RK4 through the packet, analytic continuation after it) and complex
coupling c = i pi + kappa_h / 2, at finite windows, against CF-B with the same c.
    A_TT  = [W_TT(-eta) / |R0(-eta)|] / [W_TT(+eta) / |R0(+eta)|]       (= H ratio; |G| cancels)
    A_CF  = H_CF-B(-eta) / H_CF-B(+eta)
    raw   = W_TT(-eta) / W_TT(+eta)   (what a weight measurement gives: A times |R0(-eta)| / |R0(+eta)|)
Pre-registered rule (notebook Entry 78 + amendments): (1) kappa_h = 0: A_TT = 1 within 0.5 % everywhere, else the
setup is suspect; (2) kappa_h = 0.93 kappa: |A_TT / A_CF - 1| > 0.5 % at any (f, T_w) rejects the prediction;
(3) size and dA/d(kappa_h/kappa) reported; two packet spectra (0.9 and 10 cycles) must give the same A to 0.5 %.
Usage: python e78_side_asymmetry.py [run] [fig]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e75_window_cf import packet_specs, setup, tt_weight  # noqa: E402

OUT = HERE / "e78_side_asymmetry.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e78_side_asymmetry.png"
RATIO = 0.93
F = np.array([20.0, 30.0, 40.0, 60.0, 87.0, 100.0, 150.0])
T_NS = (1.4, 14.0, 52.0, 300.0)
# First run at S/100 (peak strain 1e-5) failed the pre-registered S-convergence check (A moved 0.54 % to S/1000)
# and rule 1 by the same amount: the +-eta ratio is second-order sensitive. Rerun in the linear regime, S/1000.
S_LIN = 1000.0


def root_c(eta, kappa, kappa_h):
    c = 1j * np.pi + 0.5 * kappa_h
    r = np.roots([c, 1j * (eta + 1j * kappa), c])
    return r[np.argmin(np.abs(r))], c


def integrate_packet(st, kappa_h, per_cycle=64, s_cut=30.0):
    """As e75.integrate_packet, with eta = +-W (both signs) and complex c."""
    kappa, S, tau, q0 = st["kappa"], st["S"], st["tau"], st["q0"]
    eta = np.concatenate([st["eta"], -st["eta"]])
    ie = 1j * (eta + 1j * kappa)
    rc = [root_c(e, kappa, kappa_h) for e in eta]
    r1 = np.array([r for r, _ in rc])
    c = rc[0][1]
    r2 = 1 / r1
    s_p = s_cut * tau
    per_unit = per_cycle * max(np.abs(eta).max(), q0) / (2 * np.pi)
    n = int(np.ceil(s_p * per_unit))
    h = s_p / n
    s = np.arange(n + 1) * h
    g = lambda x: np.where(x > 0, np.exp(-np.clip(x, 0, None) / tau) * np.sin(q0 * x), 0.0)  # noqa: E731

    def f(si, r):
        e = np.exp(1j * S * g(si))
        return c * np.conj(e) + (ie + c * e * r) * r

    R = r1.copy()
    out = np.empty((n + 1, eta.size), complex)
    out[0] = R
    for i in range(n):
        k1 = f(s[i], R)
        k2 = f(s[i] + h / 2, R + h / 2 * k1)
        k3 = f(s[i] + h / 2, R + h / 2 * k2)
        k4 = f(s[i] + h, R + h * k3)
        R = R + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = R
    return dict(s=s, R=out, r1=r1, r2=r2, s_p=s_p, h=h, Lam=ie + 2 * c * r1, eta=eta)


def cf_H(st, kappa_h, T, sign):
    kappa, S = st["kappa"], st["S"]
    out = []
    for W in st["eta"]:
        eta = sign * W
        R0, c = root_c(eta, kappa, kappa_h)
        lam = 1j * (eta + 1j * kappa) + 2 * c * R0
        G, L = packet_specs(st["q0"], st["tau"], eta, lam, T)
        Ft = {Wk: 1j * c * S * (1 - R0**2) * (G(Wk) - np.exp((lam + 1j * Wk) * T) * L) / (lam + 1j * Wk) for Wk in (W, -W)}
        Iw = np.conj(R0) * Ft[W] + R0 * np.conj(Ft[-W])
        out.append(abs(Iw) / (abs(G(W)) * abs(R0) * np.pi * S / kappa))
    return np.array(out)


def one(E, kappa_h_ratio, ncyc=0.9, s_lin=S_LIN, per_cycle=64):
    st = setup(E, ncyc=ncyc, s_lin=s_lin)
    kh = kappa_h_ratio * st["kappa"]
    run = integrate_packet(st, kh, per_cycle=per_cycle)
    st2 = dict(st, eta=run["eta"])  # both signs
    nf = F.size
    rows = []
    for T_ns in T_NS:
        T = T_ns * 1e-9 / st["t_unit"]
        w, Ta = tt_weight(st2, run, T)
        r0 = np.abs(run["r1"])
        Hp, Hm = w[:nf] / r0[:nf], w[nf:] / r0[nf:]
        A_tt = Hm / Hp
        raw = w[nf:] / w[:nf]
        A_cf = cf_H(st, kh, Ta, -1) / cf_H(st, kh, Ta, +1)
        rows.append(dict(T_ns=T_ns, A_tt=A_tt.tolist(), A_cf=A_cf.tolist(), raw=raw.tolist(), r0_ratio=(r0[nf:] / r0[:nf]).tolist()))
    return rows


def section_run():
    t0 = time.time()
    out = dict(f_GHz=F.tolist(), T_ns=T_NS, ratio=RATIO, cases={})
    for E in (17.48, 8.0):
        # rule 1: symmetry without kappa_h
        rows0 = one(E, 0.0)
        w0 = max(np.max(np.abs(np.array(r["A_tt"]) - 1)) for r in rows0)
        print(f"Si 400 {E:g} keV, kappa_h = 0: max |A_TT - 1| over f, T_w = {w0 * 100:.3f} %  (rule 1: < 0.5 %)")
        # rule 2: prediction
        rows = one(E, RATIO)
        w2 = max(np.max(np.abs(np.array(r["A_tt"]) / np.array(r["A_cf"]) - 1)) for r in rows)
        print(f"Si 400 {E:g} keV, kappa_h = {RATIO} kappa: max |A_TT/A_CF - 1| = {w2 * 100:.3f} %  (rule 2: < 0.5 %)")
        for r in rows:
            print(f"   T_w = {r['T_ns']:5g} ns  A_TT  " + " ".join(f"{f:g}:{v:.4f}" for f, v in zip(F, r["A_tt"], strict=True)))
            print(f"                   A_CF  " + " ".join(f"{f:g}:{v:.4f}" for f, v in zip(F, r["A_cf"], strict=True)))
            print(f"                   raw W " + " ".join(f"{f:g}:{v:.4f}" for f, v in zip(F, r["raw"], strict=True))
                  + "   |R0| ratio " + " ".join(f"{v:.4f}" for v in r["r0_ratio"]))
        # spectrum cancellation: second packet
        rows10 = one(E, RATIO, ncyc=10.0)
        wsp = max(np.max(np.abs(np.array(a["A_tt"]) / np.array(b["A_tt"]) - 1)) for a, b in zip(rows, rows10, strict=True))
        print(f"   two packet spectra (0.9 vs 10 cycles): max |A_TT ratio - 1| = {wsp * 100:.3f} %  (< 0.5 %)")
        # sensitivity from the closed form, steady, at 20 and 100 GHz
        st = setup(E)
        T = 300e-9 / st["t_unit"]
        dA = []
        for rr in (RATIO - 0.05, RATIO + 0.05):
            kh = rr * st["kappa"]
            dA.append(cf_H(st, kh, T, -1) / cf_H(st, kh, T, +1))
        sens = (dA[1] - dA[0]) / 0.1
        print("   dA/d(kappa_h/kappa), steady: " + " ".join(f"{f:g}:{v:.3f}" for f, v in zip(F, sens, strict=True)))
        out["cases"][f"E{E:g}"] = dict(sym=w0, pred=w2, spectra=wsp, rows0=rows0, rows=rows, rows10=rows10, sens=sens.tolist())
    # convergence at 17.48 keV, 14 ns
    a = one(17.48, RATIO)[1]["A_tt"]
    b = one(17.48, RATIO, s_lin=10 * S_LIN)[1]["A_tt"]
    c_ = one(17.48, RATIO, per_cycle=128)[1]["A_tt"]
    print(f"convergence (17.48 keV, 14 ns): S/{S_LIN:g} -> S/{10 * S_LIN:g} max |dA| = {np.max(np.abs(np.array(a) - np.array(b))) * 100:.3f} %, step 64 -> 128 {np.max(np.abs(np.array(a) - np.array(c_))) * 100:.3f} %")
    out["conv"] = dict(S=float(np.max(np.abs(np.array(a) - np.array(b)))), step=float(np.max(np.abs(np.array(a) - np.array(c_)))))
    ok1 = all(out["cases"][k]["sym"] < 0.005 for k in out["cases"])
    ok2 = all(out["cases"][k]["pred"] < 0.005 for k in out["cases"])
    out["verdict"] = ("setup suspect" if not ok1 else ("prediction confirmed in the model" if ok2 else "prediction rejected"))
    print("verdict:", out["verdict"], f"({time.time() - t0:.0f} s)")
    json.dump(out, open(OUT, "w"), indent=1)


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    f = np.array(d["f_GHz"])
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True, sharey=True)
    for ax, E in zip(axs, ("17.48", "8"), strict=True):
        c = d["cases"][f"E{E}"]
        for r, col in zip(c["rows"], ("C0", "C1", "C3", "k"), strict=True):
            ax.plot(f, 100 * (np.array(r["A_cf"]) - 1), "-", color=col, lw=1.4, label=f"closed form, {r['T_ns']:g} ns")
            ax.plot(f, 100 * (np.array(r["A_tt"]) - 1), "o", color=col, ms=4, mfc="none")
        ax.plot(f, 100 * (np.array(c["rows"][-1]["raw"]) - 1), "--", color="0.5", lw=1, label="raw weight ratio, 300 ns (incl. |R0| ratio)")
        ax.axhline(0, color="0.6", lw=0.6)
        ax.set(xscale="log", xlabel="phonon frequency (GHz)", title=f"Si 400, {E} keV, $\\chi_{{ih}} = 0.93\\,\\chi_{{i0}}$")
        ax.grid(color="#e4e3df", lw=0.6)
    axs[0].set_ylabel("$H(-\\eta)/H(+\\eta) - 1$ (%)  [lines: closed form, circles: TT]")
    axs[0].legend(fontsize=7, frameon=False)
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
