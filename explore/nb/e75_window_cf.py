"""Entry 75 (T2a): the crystal factor H as a closed form in the window length T_w, tested against TT.

Entry 71 showed that the Entry 65 crystal factor H(f) is the infinite-window limit and gave a hand estimate
for a 14 ns window. Here the first-order (Entry 59) response of a packet that enters the surface at t = 0 is
written with the window kept:
    delta(t) = -i c (1 - R0^2) S int_0^t e^{Lambda (t - s)} g(s) ds,     Lambda = i (eta + i kappa) + 2 c R0,
    F_T(W)   = int_0^T delta e^{iWt} dt = i c (1 - R0^2) S [G_T(W) - e^{(Lambda + iW) T} L(Lambda)] / (Lambda + iW),
    G_T(W)   = int_0^T g e^{iWs} ds,   L(Lambda) = int_0^T g e^{-Lambda s} ds,
    I~(eta)  = conj(R0) F_T(eta) + R0 conj(F_T(-eta))   (Fourier component of |R|^2 - |R0|^2 at eta).
CF-B keeps G and L exactly (analytic for the damped-sinusoid packet of Entry 58/66).
CF-A is the packet-free crystal factor: a packet short against min(1/kappa, x/pi) and inside the window has
L = G(-eta), G_T = G, which leaves the Entry 65 cos-mode form with each kernel 1/(Lambda - i s W) replaced by
[1 - e^{(Lambda - i s W) T}] / (Lambda - i s W). Its leading term is H(f; T) = H_inf(f) |1 - e^{-(kappa + i pi/x) T}|.

TT: the Entry 66 semi-infinite trick (one Riccati ODE in s = t - zeta per eta). The packet is integrated by RK4
up to s_p where g < 1e-13; beyond s_p the crystal is unstrained and the Riccati has constant coefficients, so
u = (R - r1)/(R - r2) = u_p e^{Lambda (s - s_p)} exactly (r1 = R0, r2 = 1/R0 the two roots). |R|^2 - |R0|^2 is
then a sum of exponentials whose Fourier integral over [s_p, T] is analytic: every window costs nothing beyond
the one RK4 through the packet. The continuation is checked against a plain RK4 run at 14 ns (Entry 71).

Weight: |Fourier component at |eta| of |R(t)|^2 - |R0|^2 over [0, T]|, divided by |G(eta)| |R0| pi alpha S / kappa
(the kinematic thin-absorber law, Entry 59/65), which makes TT, CF-A and CF-B directly comparable as H.
The window-mean-subtracted weight of Entry 71 is also computed for the 140 ps row (report only).

Pre-registered rule (notebook Entry 75): |H_TT / H_CF-A - 1| > 1 % at any (f, T_w) for the Entry 66 packet
(convergence checks passing) rejects CF-A; CF-B is then tested with the same rule.
Usage: python e75_window_cf.py [run] [fig]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e65_si400 import params, root  # noqa: E402

V_L = 8433.0
F = np.array([20.0, 30.0, 40.0, 60.0, 87.0, 100.0, 150.0])
T_NS = np.array([0.14, 0.5, 1.4, 3.0, 5.0, 8.0, 14.0, 22.0, 30.0, 50.0, 80.0, 140.0, 300.0])
OUT = HERE / "e75_window_cf.json"


def setup(E, f0=87.0, ncyc=0.9, s_lin=100.0):
    p = params(E)
    t_unit = p["Lambda_um"] * 1e-6 / V_L
    g_lambda = 2 * np.pi / 1.3577e-10 * p["Lambda_um"] * 1e-6
    q0 = 2 * np.pi * f0 * 1e9 * t_unit
    eta = 2 * np.pi * F * 1e9 * t_unit
    kappa = 2 * np.pi * p["g"]
    S = 1e-3 * g_lambda / q0 / s_lin  # peak strain 1e-3 / s_lin
    tau = ncyc * 2 * np.pi / q0
    return dict(p=p, t_unit=t_unit, q0=q0, eta=eta, kappa=kappa, S=S, tau=tau, ncyc=ncyc, E=E)


def packet_specs(q0, tau, eta, Lam, T=np.inf):
    """G_T(W) = int_0^T e^{-s/tau} sin(q0 s) e^{iWs} ds and L_T(Lambda) = int_0^T ... e^{-Lambda s} ds, analytic."""

    def lap(p):
        if np.isinf(T):
            return q0 / (p**2 + q0**2)
        return (q0 - np.exp(-p * T) * (q0 * np.cos(q0 * T) + p * np.sin(q0 * T))) / (p**2 + q0**2)

    G = lambda W: lap(1 / tau - 1j * W)  # noqa: E731
    return G, lap(1 / tau + Lam)


def closed_forms(st, T):
    """H_CF-A and H_CF-B (cos-mode and packet-exact first order), window T (units Lambda/v_L)."""
    eta, kappa, S = st["eta"], st["kappa"], st["S"]
    c = 1j * np.pi
    A, B = [], []
    for e in eta:
        R0 = root(1.0, e, kappa)
        Lam = 1j * (e + 1j * kappa) + 2 * c * R0
        norm = abs(R0) * np.pi * S / kappa
        # CF-A: Entry 65 cos-mode form, kernels truncated at T
        D = {s: 1j * c * S * (1 - R0**2) / 2 * (1 - np.exp((Lam - 1j * s * e) * T)) / (Lam - 1j * s * e) for s in (1, -1)}
        A.append(abs(np.conj(R0) * D[-1] + R0 * np.conj(D[1])) / (norm / 2))
        # CF-B: packet kept, spectra truncated at the window (G_T at +-eta, L_T at Lambda)
        G, L = packet_specs(st["q0"], st["tau"], e, Lam, T)
        Ft = {}
        for W in (e, -e):
            Ft[W] = 1j * c * S * (1 - R0**2) * (G(W) - np.exp((Lam + 1j * W) * T) * L) / (Lam + 1j * W)
        I = np.conj(R0) * Ft[e] + R0 * np.conj(Ft[-e])
        B.append(abs(I) / (abs(G(e)) * norm))  # same |G_T(eta)| normalisation as H_TT
    return np.array(A), np.array(B)


def integrate_packet(st, per_cycle=64, s_cut=30.0):
    """RK4 through the packet, s in [0, s_p], s_p = s_cut tau. Returns s grid, R(s, eta), and roots."""
    eta, kappa, S = st["eta"], st["kappa"], st["S"]
    tau, q0 = st["tau"], st["q0"]
    c = 1j * np.pi
    ie = 1j * (eta + 1j * kappa)
    r1 = np.array([root(1.0, e, kappa) for e in eta])
    r2 = 1 / r1
    s_p = s_cut * tau
    per_unit = per_cycle * max(eta.max(), q0) / (2 * np.pi)
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
    return dict(s=s, R=out, r1=r1, r2=r2, s_p=s_p, h=h, Lam=ie + 2 * c * r1)


def tt_weight(st, run, T, omega=None, mean_sub=False):
    """|int_0^T (|R|^2 - |R0|^2) e^{i omega t} dt| for omega = eta (default). Returns (weight, T_actual)."""
    s, R, r1, r2, Lam = run["s"], run["R"], run["r1"], run["r2"], run["Lam"]
    eta = st["eta"] if omega is None else omega
    I = np.abs(R) ** 2 - np.abs(r1) ** 2
    if T <= run["s_p"]:
        n = int(round(T / run["h"]))
        T = n * run["h"]
        w = np.full(n + 1, run["h"])
        w[[0, -1]] /= 2
        tot = np.sum(w[:, None] * I[: n + 1] * np.exp(1j * eta[None, :] * s[: n + 1, None]), axis=0)
        if mean_sub:
            tot -= np.sum(w[:, None] * I[: n + 1], axis=0) / T * (np.exp(1j * eta * T) - 1) / (1j * eta)
        return np.abs(tot), T
    w = np.full(s.size, run["h"])
    w[[0, -1]] /= 2
    tot = np.sum(w[:, None] * I * np.exp(1j * eta[None, :] * s[:, None]), axis=0)
    mean_part = np.sum(w[:, None] * I, axis=0)
    # analytic continuation: delta(s) = (r1 - r2) sum_k u_p^k e^{k Lambda (s - s_p)}
    u_p = (R[-1] - r1) / (R[-1] - r2)
    assert np.all(np.abs(u_p) < 0.5), np.abs(u_p).max()
    K = int(np.ceil(np.log(1e-15) / np.log(np.abs(u_p).max()))) + 1
    k = np.arange(1, K + 1)
    d = (r1 - r2)[None, :] * u_p[None, :] ** k[:, None]  # (K, n_eta)
    sp, dT = run["s_p"], T - run["s_p"]

    def fint(a, om):
        """int_{s_p}^T e^{a (s - s_p)} e^{i om s} ds."""
        b = a + 1j * om
        return np.exp(1j * om * sp) * (np.exp(b * dT) - 1) / b

    lin = np.zeros(eta.size, complex)
    mean_lin = np.zeros(eta.size, complex)
    for j in range(eta.size):
        a = k * Lam[j]
        lin[j] = np.conj(r1[j]) * np.sum(d[:, j] * fint(a, eta[j])) + r1[j] * np.sum(np.conj(d[:, j]) * fint(np.conj(a), eta[j]))
        mean_lin[j] = np.conj(r1[j]) * np.sum(d[:, j] * fint(a, 0.0)) + r1[j] * np.sum(np.conj(d[:, j]) * fint(np.conj(a), 0.0))
        # |delta|^2 term, second order (kept so that the continuation is exact, not linearised)
        aa = (k[:, None] * Lam[j] + k[None, :] * np.conj(Lam[j])).ravel()
        dd = (d[:, j][:, None] * np.conj(d[:, j])[None, :]).ravel()
        lin[j] += np.sum(dd * fint(aa, eta[j]))
        mean_lin[j] += np.sum(dd * fint(aa, 0.0))
    tot = tot + lin
    if mean_sub:
        tot -= (mean_part + mean_lin) / T * (np.exp(1j * eta * T) - 1) / (1j * eta)
    return np.abs(tot), T


def h_norm(st, T=np.inf):
    """Kinematic normalisation |G_T(eta)| |R0| pi S / kappa, per eta (strain spectrum over the window)."""
    eta, kappa, S = st["eta"], st["kappa"], st["S"]
    out = []
    for e in eta:
        R0 = root(1.0, e, kappa)
        G, _ = packet_specs(st["q0"], st["tau"], e, 0.0, T)
        out.append(abs(G(e)) * abs(R0) * np.pi * S / kappa)
    return np.array(out)


def plain_rk4(st, T, per_cycle=64):
    """Entry 66/71 style: RK4 all the way to T (check of the analytic continuation)."""
    from e66_harmonic import surface_R

    eta, kappa, S = st["eta"], st["kappa"], st["S"]
    tau, q0 = st["tau"], st["q0"]
    g = lambda x: np.where(x > 0, np.exp(-np.clip(x, 0, None) / tau) * np.sin(q0 * x), 0.0)  # noqa: E731
    per_unit = per_cycle * max(eta.max(), q0) / (2 * np.pi)
    t, R = surface_R(1.0, kappa, eta, lambda x: S * g(x), T, per_unit)
    I = np.abs(R) ** 2 - np.abs(R[0]) ** 2
    w = np.full(t.size, t[1] - t[0])
    w[[0, -1]] /= 2
    return np.abs(np.sum(w[:, None] * I * np.exp(1j * eta[None, :] * t[:, None]), axis=0))


def section_run():
    out = dict(f_GHz=F.tolist(), T_ns=T_NS.tolist(), cases={})
    t0 = time.time()
    # --- checks on the main case ---
    st = setup(17.48)
    run = integrate_packet(st)
    print(f"Si 400 17.48 keV: Lambda = {st['p']['Lambda_um']:.2f} um, t_unit = {st['t_unit'] * 1e9:.3f} ns, "
          f"kappa = {st['kappa']:.4f}, S = {st['S']:.2e}, s_p = {run['s_p'] * st['t_unit'] * 1e12:.0f} ps, "
          f"RK4 steps = {run['s'].size - 1} ({time.time() - t0:.1f} s)")
    T14 = 14e-9 / st["t_unit"]
    w_cont, _ = tt_weight(st, run, T14)
    w_plain = plain_rk4(st, T14)
    print("check continuation vs plain RK4 at 14 ns: ratio = " + " ".join(f"{a / b:.5f}" for a, b in zip(w_cont, w_plain, strict=True)))
    out["check_continuation_14ns"] = (w_cont / w_plain).tolist()
    w_ms, _ = tt_weight(st, run, T14, mean_sub=True)
    hn = h_norm(st)
    r = (w_ms / hn) / (w_ms / hn)[list(F).index(100.0)]
    print("check Entry 71 (mean-subtracted, /H(100)) at 14 ns: " + " ".join(f"{f:g}:{v:.3f}" for f, v in zip(F, r, strict=True)))
    out["check_e71_14ns"] = r.tolist()
    run128 = integrate_packet(st, per_cycle=128)
    w128, _ = tt_weight(st, run128, T14)
    print("check step 64 -> 128 at 14 ns: ratio = " + " ".join(f"{a / b:.5f}" for a, b in zip(w_cont, w128, strict=True)))
    out["check_step_14ns"] = (w_cont / w128).tolist()
    st1k = setup(17.48, s_lin=1000.0)
    run1k = integrate_packet(st1k)
    w1k, _ = tt_weight(st1k, run1k, T14)
    print("check S/100 -> S/1000 at 14 ns: ratio of H = " + " ".join(f"{a / b:.5f}" for a, b in zip(w_cont / hn, w1k / h_norm(st1k), strict=True)))
    out["check_S_14ns"] = ((w_cont / hn) / (w1k / h_norm(st1k))).tolist()
    # --- scan: energies x packets x windows ---
    worst = {}
    for E in (17.48, 8.0):
        for ncyc in (0.9, 10.0, 100.0):
            st = setup(E, ncyc=ncyc)
            t1 = time.time()
            run = integrate_packet(st)
            hn = h_norm(st)
            key = f"E{E:g}_ncyc{ncyc:g}"
            rows = []
            wa = wb = 0.0
            for T_ns in T_NS:
                T = T_ns * 1e-9 / st["t_unit"]
                w, Ta = tt_weight(st, run, T)
                A, B = closed_forms(st, Ta)
                Htt = w / h_norm(st, Ta)
                rows.append(dict(T_ns=T_ns, T=Ta, H_tt=Htt.tolist(), H_A=A.tolist(), H_B=B.tolist()))
                wa = max(wa, np.abs(Htt / A - 1).max())
                wb = max(wb, np.abs(Htt / B - 1).max())
            worst[key] = (wa, wb)
            out["cases"][key] = dict(E=E, ncyc=ncyc, tau_ns=st["tau"] * st["t_unit"] * 1e9, kappa=st["kappa"], S=st["S"], rows=rows)
            print(f"{key}: tau = {st['tau'] * st['t_unit'] * 1e12:.0f} ps, steps = {run['s'].size - 1}, "
                  f"worst |H_TT/H_CF-A - 1| = {wa * 100:.2f} %, |H_TT/H_CF-B - 1| = {wb * 100:.3f} % ({time.time() - t1:.1f} s)")
            if ncyc == 0.9:
                i100 = list(F).index(100.0)
                print("   T_ns     " + "".join(f"{f:>8g}" for f in F) + "   | H(20)/H(100): TT  CF-A  steady")
                for r_ in rows:
                    Htt, A = np.array(r_["H_tt"]), np.array(r_["H_A"])
                    print(f"   {r_['T_ns']:6g}  TT " + "".join(f"{v:8.4f}" for v in Htt)
                          + f"   | {Htt[0] / Htt[i100]:.3f} {A[0] / A[i100]:.3f}")
                    print("          CF-A " + "".join(f"{v:8.4f}" for v in A)
                          + "   ratio " + "".join(f"{v:6.3f}" for v in Htt / A))
    out["worst"] = worst
    # --- boundaries from the closed form that passes the rule (CF-B, 17.48 keV), then confirmed by TT ---
    st = setup(17.48)
    run = integrate_packet(st)
    Tg = np.geomspace(0.1, 400, 2000) * 1e-9 / st["t_unit"]
    i20, i100 = 0, list(F).index(100.0)
    ratio = np.array([(lambda B: B[i20] / B[i100])(closed_forms(st, T)[1]) for T in Tg])
    steady = closed_forms(st, 1e4)[1]
    r_inf = steady[i20] / steady[i100]
    bounds = {}
    for name, target in (("T_flat", 0.9), ("T_half", 0.5 * (1 + r_inf)), ("T_90", r_inf + 0.1 * (1 - r_inf))):
        j = np.argmax(ratio < target)
        T = Tg[j]
        w, Ta = tt_weight(st, run, T)
        _, B = closed_forms(st, Ta)
        Htt = w / h_norm(st, Ta)
        bounds[name] = dict(T_ns=T * st["t_unit"] * 1e9, target=target, cf=B[i20] / B[i100], tt=Htt[i20] / Htt[i100])
        print(f"{name}: ratio target {target:.3f} at T = {T * st['t_unit'] * 1e9:.1f} ns  (CF-B {B[i20] / B[i100]:.3f}, TT {Htt[i20] / Htt[i100]:.3f})")
    out["bounds"] = bounds
    out["ratio_curve"] = dict(T_ns=(Tg * st["t_unit"] * 1e9).tolist(), ratio=ratio.tolist(), steady=r_inf)
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT.name} ({time.time() - t0:.0f} s)")


def section_posthoc():
    """Written after the run (not pre-registered). (1) Where CF-A fails, is it the packet's group delay
    t_g(eta) = d arg G / dW, i.e. the window effectively starting t_g after the front enters? Test
    CF-A with T -> T - t_g (first-order Taylor of CF-B in the packet duration). (2) Is the 0.8 % of
    CF-B at ncyc = 100 second order in S? Rerun at S/1000."""
    res = {}
    st = setup(17.48)
    run = integrate_packet(st)
    rows = []
    for T_ns in T_NS[:5]:
        T = T_ns * 1e-9 / st["t_unit"]
        w, Ta = tt_weight(st, run, T)
        Htt = w / h_norm(st, Ta)
        A = []
        for e in st["eta"]:
            G, _ = packet_specs(st["q0"], st["tau"], e, 0.0)
            dW = 1e-4 * e
            tg = (np.angle(G(e + dW)) - np.angle(G(e - dW))) / (2 * dW)
            st1 = dict(st, eta=np.array([e]))
            A.append(closed_forms(st1, Ta - tg)[0][0])
        A = np.array(A)
        rows.append(dict(T_ns=T_ns, ratio=(Htt / A).tolist()))
        print(f"post-hoc CF-A with T - t_g, {T_ns:5g} ns: TT/CF-A' = " + " ".join(f"{v:.4f}" for v in Htt / A))
    res["group_delay"] = rows
    for ncyc in (10.0, 100.0):
        wb = {}
        for s_lin in (100.0, 1000.0):
            st = setup(17.48, ncyc=ncyc, s_lin=s_lin)
            run = integrate_packet(st)
            worst = 0.0
            for T_ns in T_NS:
                T = T_ns * 1e-9 / st["t_unit"]
                w, Ta = tt_weight(st, run, T)
                _, B = closed_forms(st, Ta)
                worst = max(worst, np.abs(w / h_norm(st, Ta) / B - 1).max())
            wb[s_lin] = worst
        print(f"post-hoc ncyc = {ncyc:g}: worst |H_TT/H_CF-B - 1| at S/100 = {wb[100.0] * 100:.3f} %, at S/1000 = {wb[1000.0] * 100:.3f} %")
        res[f"cfb_ncyc{ncyc:g}"] = wb
    d = json.load(open(OUT))
    d["posthoc"] = res
    json.dump(d, open(OUT, "w"), indent=1)


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    c = d["cases"]["E17.48_ncyc0.9"]
    f = np.array(d["f_GHz"])
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    ax = axs[0]
    Tn = np.array([r["T_ns"] for r in c["rows"]])
    for j, fj in enumerate(f):
        if fj not in (20, 40, 100, 150):
            continue
        Htt = np.array([r["H_tt"][j] for r in c["rows"]])
        HA = np.array([r["H_A"][j] for r in c["rows"]])
        HB = np.array([r["H_B"][j] for r in c["rows"]])
        (ln,) = ax.plot(Tn, HB, "-", lw=1.5, label=f"{fj:g} GHz")
        ax.plot(Tn, Htt, "o", color=ln.get_color(), ms=4, mfc="none")
    ax.set(xscale="log", xlabel="window $T_w$ (ns)", ylabel="crystal factor $H(f; T_w)$", title="Si 400, 17.48 keV: first-order closed form (lines), TT (circles)")
    ax.axvline(22, color="0.6", lw=0.8, ls=":")
    ax.text(23, 0.05, r"$1/\kappa$", color="0.4", fontsize=9)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(color="#e4e3df", lw=0.6)
    ax = axs[1]
    rc = d["ratio_curve"]
    ax.plot(rc["T_ns"], rc["ratio"], "-", color="C3", lw=1.5, label="closed form (CF-B)")
    i100 = list(f).index(100.0)
    ax.plot(Tn, [r["H_tt"][0] / r["H_tt"][i100] for r in c["rows"]], "o", color="C3", ms=4, mfc="none", label="TT")
    ax.axhline(rc["steady"], color="0.5", lw=0.8, ls="--")
    ax.text(150, rc["steady"] + 0.012, "steady (Entry 65)", color="0.4", fontsize=9, ha="center")
    for k, (name, b) in enumerate(d["bounds"].items()):
        ax.axvline(b["T_ns"], color="0.7", lw=0.8, ls=":")
        ax.text(b["T_ns"], 1.03 - 0.045 * k, f"{name.replace('T_', '')} {b['T_ns']:.0f} ns", fontsize=8, color="0.4", ha="right")
    ax.set(xscale="log", xlabel="window $T_w$ (ns)", ylabel="$H(20\\,\\mathrm{GHz}; T_w)\\,/\\,H(100\\,\\mathrm{GHz}; T_w)$", ylim=(0.45, 1.06), title="where the Entry 65 factor 2 comes back")
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    ax.grid(color="#e4e3df", lw=0.6)
    fig.savefig(HERE.parents[1] / "figures" / "paper" / "e75_window_H.png", dpi=200)
    print("wrote figures/paper/e75_window_H.png")


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run"]:
        globals()[f"section_{sec}"]()
