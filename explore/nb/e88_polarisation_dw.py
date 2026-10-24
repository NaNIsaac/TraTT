"""Entry 88 (T2r): polarisation (pi component, C = |cos 2 theta_B|) and the static Debye-Waller factor as rescalings.

pi component: both couplings c -> C c, kappa unchanged: alpha_pi = C alpha, g_pi = g / C, x_pi = x / C, nu_pi = C nu.
Unpolarised beam: F_mix = F_sigma + F_pi (deterministic, same strain), H_mix = |F_sigma + F_pi| / (N_sigma + N_pi).
Checks (pre-registered, notebook Entry 88): (i) pure pi = sigma closed form with (g/C, x/C) against the TT run at
c = i pi C (identity, one point); (ii) H_mix closed-form sum vs TT sum (sigma run + pi run) to 0.3 % at every (f, T_w)
for Si 440 / 17.48 keV Bragg (longitudinal [110] wave) and Si 400 / 8 keV. Reported: H_pi/H_sigma, H_mix/H_sigma,
window numbers for sigma, pi, unpolarised; nu_pi/nu_sigma = C; e^{-M} values.
Usage: python e88_polarisation_dw.py
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e65_si400 import SI400, params, root  # noqa: E402
from e75_window_cf import packet_specs, setup, tt_weight  # noqa: E402
from e81_asym_map import integrate  # noqa: E402

OUT = HERE / "e88_polarisation_dw.json"
V_L110, V_L100 = 9136.0, 8433.0
HC = 12.398419843
# Si 440 at 17.48 keV (X0h, e20): theta_B, |chi_rh|, chi_i0 ; Bragg extinction length lambda sin(theta)/|chi_rh|
SI440 = (21.679, 0.12130e-5, 0.16059e-7)


def cf_complex_c(st, c, eta, W, T):
    """CF-B (Entry 75) with an explicit coupling c, at detuning eta, Fourier frequency W (both signs summed), window T."""
    kappa, S = st["kappa"], st["S"]
    r = np.roots([c, 1j * (eta + 1j * kappa), c])
    R0 = r[np.argmin(np.abs(r))]
    lam = 1j * (eta + 1j * kappa) + 2 * c * R0
    G, L = packet_specs(st["q0"], st["tau"], eta, lam, T)
    Ft = {Wk: 1j * c * S * (1 - R0**2) * (G(Wk) - np.exp((lam + 1j * Wk) * T) * L) / (lam + 1j * Wk) for Wk in (W, -W)}
    I = np.conj(R0) * Ft[W] + R0 * np.conj(Ft[-W])
    N = abs(G(W)) * abs(R0) * abs(c) * S / kappa  # kinematic law with |c| = pi alpha C
    return I, N


def fourier_complex(st, run, T):
    """Complex Fourier components at omega = eta (per eta), e82 formula, with the e81 run dict."""
    from e82_divergence import fourier_complex as fc

    return fc(st, run, T, st["eta"])


def case(label, E, theta_B, xrh, xi0, lam_um, v, f_GHz, T_ns_list, S_lin=1000.0):
    th = np.radians(theta_B)
    C = abs(np.cos(2 * th))
    lam = lam_um * 1e-6
    t_unit = lam / v
    g = xi0 / xrh
    base = setup(E, s_lin=S_lin) if E in SI400 else None
    st = dict(t_unit=t_unit, kappa=2 * np.pi * g, q0=2 * np.pi * 87e9 * t_unit, tau=0.9 * 2 * np.pi / (2 * np.pi * 87e9 * t_unit))
    st["S"] = base["S"] if base else 7.14e-6  # peak strain ~1e-6 scale irrelevant at first order (ratios only)
    st["eta"] = 2 * np.pi * np.array(f_GHz) * 1e9 * t_unit
    res = dict(label=label, C=C, g=g, f_GHz=list(f_GHz), x=(st["eta"] / (2 * np.pi)).tolist(), rows=[], worst=0.0)
    c_s, c_p = 1j * np.pi, 1j * np.pi * C
    run_s = integrate(st, c_s, c_s)
    run_p = integrate(st, c_p, c_p)
    for T_ns in T_ns_list:
        T = T_ns * 1e-9 / t_unit
        Fs, Ta = fourier_complex(st, run_s, T)
        Fp, _ = fourier_complex(st, run_p, T)
        row = dict(T_ns=T_ns, H_s=[], H_p=[], H_mix=[], Hp_over_Hs=[], Hmix_over_Hs=[], tt_over_cf_mix=[], tt_over_cf_pi=[])
        for j, e in enumerate(st["eta"]):
            Is, Ns = cf_complex_c(st, c_s, e, e, Ta)
            Ip, Np = cf_complex_c(st, c_p, e, e, Ta)
            Hs, Hp, Hm = abs(Is) / Ns, abs(Ip) / Np, abs(Is + Ip) / (Ns + Np)
            Hm_tt = abs(Fs[j] + Fp[j]) / (Ns + Np)
            Hp_tt = abs(Fp[j]) / Np
            row["H_s"].append(Hs)
            row["H_p"].append(Hp)
            row["H_mix"].append(Hm)
            row["Hp_over_Hs"].append(Hp / Hs)
            row["Hmix_over_Hs"].append(Hm / Hs)
            row["tt_over_cf_mix"].append(Hm_tt / Hm)
            row["tt_over_cf_pi"].append(Hp_tt / Hp)
            res["worst"] = max(res["worst"], abs(Hm_tt / Hm - 1))
        res["rows"].append(row)
        print(f"  {label} T_w = {T_ns:g} ns: " + " | ".join(
            f"{f:g} GHz: H_s {row['H_s'][j]:.3f}, H_pi {row['H_p'][j]:.3f} ({row['Hp_over_Hs'][j]:.3f} x), H_mix {row['H_mix'][j]:.3f} ({row['Hmix_over_Hs'][j]:.3f} x), TT/CF mix {row['tt_over_cf_mix'][j]:.4f}, pi {row['tt_over_cf_pi'][j]:.4f}"
            for j, f in enumerate(f_GHz)), flush=True)
    # window numbers for sigma, pi, unpolarised: H(f_lo)/H(f_hi) against T_w, Entry 75 definitions (H = |I| / N)
    Tg = np.geomspace(0.3, 600, 1200) * 1e-9 / t_unit
    e_lo, e_hi = st["eta"][0], st["eta"][-1]
    win = {}
    for pol, cc in (("sigma", c_s), ("pi", c_p), ("unpol", None)):
        def ratio_at(T):
            if cc is not None:
                a, na = cf_complex_c(st, cc, e_lo, e_lo, T)
                b, nb = cf_complex_c(st, cc, e_hi, e_hi, T)
                return (abs(a) / na) / (abs(b) / nb)
            a1, n1 = cf_complex_c(st, c_s, e_lo, e_lo, T)
            a2, n2 = cf_complex_c(st, c_p, e_lo, e_lo, T)
            b1, m1 = cf_complex_c(st, c_s, e_hi, e_hi, T)
            b2, m2 = cf_complex_c(st, c_p, e_hi, e_hi, T)
            return (abs(a1 + a2) / (n1 + n2)) / (abs(b1 + b2) / (m1 + m2))
        r_inf = ratio_at(1e4)
        rr = np.array([ratio_at(T) for T in Tg])
        w = {}
        for name, target in (("T_flat", 1 - 0.1 * (1 - r_inf)), ("T_half", 0.5 * (1 + r_inf)), ("T_90", r_inf + 0.1 * (1 - r_inf))):
            j = int(np.argmax((rr - target) * np.sign(1 - r_inf) < 0))
            w[name] = float(Tg[j] * t_unit * 1e9)
        w["steady_ratio"] = float(r_inf)
        win[pol] = w
        print(f"  {label} window law ({pol}): H({f_GHz[0]:g})/H({f_GHz[-1]:g}) steady {r_inf:.3f}; T_flat {w['T_flat']:.1f} / T_half {w['T_half']:.1f} / T_90 {w['T_90']:.1f} ns  (targets 10 % / 50 % / 90 % of the way from 1 to steady)")
    res["window"] = win
    return res


def main():
    out = dict(cases=[], dw={})
    # Si 440 / 17.48 keV Bragg, longitudinal [110] (Entry 85 addendum)
    th, xrh, xi0 = SI440
    lam440 = HC / 17.48 * 1e-10 * np.sin(np.radians(th)) / xrh
    print(f"Si 440 / 17.48 keV Bragg, longitudinal: C = {abs(np.cos(np.radians(2 * th))):.3f}, Lambda = {lam440 * 1e6:.1f} um, g = {xi0 / xrh:.4f}")
    out["cases"].append(case("Si440_17.48keV", 17.48, th, xrh, xi0, lam440 * 1e6, V_L110, (1.0, 2.0), (7.0, 300.0)))
    # Si 400 / 8 keV Bragg, longitudinal [100]
    p8 = params(8.0)
    th8, xrh8, xi08 = SI400[8.0]
    print(f"Si 400 / 8 keV Bragg: C = {abs(np.cos(np.radians(2 * th8))):.3f}, Lambda = {p8['Lambda_um']:.2f} um, g = {p8['g']:.4f}")
    out["cases"].append(case("Si400_8keV", 8.0, th8, xrh8, xi08, p8["Lambda_um"], V_L100, (20.0, 40.0, 100.0), (14.0, 300.0)))
    worst = max(c["worst"] for c in out["cases"])
    out["verdict_pol"] = "absorbed as a factor (pure pi = rescaling; unpolarised = closed-form sum, TT to 0.3 %)" if worst < 0.003 else "numbers in the limitations section"
    print(f"TT sum vs closed-form sum, worst: {worst * 100:.3f} % -> {out['verdict_pol']}")
    # Debye-Waller: e^{-M} with M = B (sin theta / lambda)^2, B = 0.47 A^2 (Si, 300 K)
    B = 0.47
    for lab, d_A in (("Si400", 1.3577), ("Si440", 0.9601), ("Si220", 1.9201)):
        M = B * (1 / (2 * d_A)) ** 2
        out["dw"][lab] = float(np.exp(-M))
    print("Debye-Waller e^{-M} (B = 0.47 A^2, 300 K): " + ", ".join(f"{k} {v:.3f}" for k, v in out["dw"].items()) + "  (already in the X0h chi_h used throughout)")
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT.name}")


if __name__ == "__main__":
    main()
