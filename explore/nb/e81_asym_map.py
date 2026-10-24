"""Entry 81 (T2g): asymmetric Bragg geometry maps exactly onto the symmetric solver; the geometry factor of H.

Two-beam Bragg with gamma_0 != |gamma_h| gives, for R = D_h / D_0 in depth units of L,
    dR/ds = c_h e^{-i phi} + i eta_c R + c_0 e^{+i phi} R^2      (retarded-time form of Entry 66, s = t - zeta),
    c_h = b_c L / |gamma_h|,  c_0 = b_c L / gamma_0,  kappa = mu L (1/gamma_0 + 1/|gamma_h|) / 2.
With R~ = R sqrt(c_0 / c_h) = R sqrt(|gamma_h| / gamma_0) this is the symmetric equation with c_eff = sqrt(c_0 c_h):
alpha_asym = L |chi_rh| / (lambda sqrt(gamma_0 |gamma_h|)), |R|^2 = |R~|^2 gamma_0 / |gamma_h|, and the crystal factor
keeps its form with g -> g (sqrt(b) + 1/sqrt(b)) / 2, b = gamma_0 / |gamma_h|.

Pre-registered rule (notebook Entry 81): semi-infinite TT with the two couplings kept distinct (Entry 66 packet at
S/100, Si 400 / 17.48 keV constants, b = 0.5, 2, 4) against the symmetric integration at c_eff with R~ rescaled, on
the same RK4 grid: if max |R_asym - sqrt(b) R~| > 1e-6 the mapping is wrong (solver extension needed). Otherwise the
geometry factor is one number, g (sqrt(b) + 1/sqrt(b)) / 2, and H_TT(b) must match the closed form at g_asym to 0.3 %
(boss's T2g criterion). Reported: H(20)/H(100), T_flat / T_half / T_90 for b = 0.5, 1, 2, 4, 10.
Usage: python e81_asym_map.py [run]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e65_si400 import SI400  # noqa: E402
from e75_window_cf import closed_forms, h_norm, setup, tt_weight  # noqa: E402

OUT = HERE / "e81_asym_map.json"
V_L = 8433.0
F = np.array([20.0, 30.0, 40.0, 60.0, 87.0, 100.0, 150.0])


def geometry(E, b):
    """Asymmetry angle chi with b = sin(theta_B + chi) / sin(theta_B - chi); Lambda_asym / Lambda_sym = sqrt(g0 gh) / sin theta_B."""
    th = np.radians(SI400[E][0])
    chi = brentq(lambda x: np.sin(th + x) / np.sin(th - x) - b, -th + 1e-6, th - 1e-6)
    g0, gh = np.sin(th + chi), np.sin(th - chi)
    return dict(chi_deg=np.degrees(chi), gamma0=g0, gammah=gh, lam_ratio=np.sqrt(g0 * gh) / np.sin(th), gfac=(np.sqrt(b) + 1 / np.sqrt(b)) / 2)


def integrate(st, c_h, c_0, per_cycle=64, s_cut=30.0):
    """RK4 through the packet with distinct couplings; returns the e75 run dict (continuation uses r1, r2, Lam)."""
    eta, kappa, S, tau, q0 = st["eta"], st["kappa"], st["S"], st["tau"], st["q0"]
    ie = 1j * (eta + 1j * kappa)
    roots = [np.roots([c_0, 1j * (e + 1j * kappa), c_h]) for e in eta]
    r1 = np.array([r[np.argmin(np.abs(r))] for r in roots])
    r2 = np.array([r[np.argmax(np.abs(r))] for r in roots])
    s_p = s_cut * tau
    per_unit = per_cycle * max(np.abs(eta).max(), q0) / (2 * np.pi)
    n = int(np.ceil(s_p * per_unit))
    h = s_p / n
    s = np.arange(n + 1) * h
    g = lambda x: np.where(x > 0, np.exp(-np.clip(x, 0, None) / tau) * np.sin(q0 * x), 0.0)  # noqa: E731

    def f(si, r):
        e = np.exp(1j * S * g(si))
        return c_h * np.conj(e) + (ie + c_0 * e * r) * r

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
    return dict(s=s, R=out, r1=r1, r2=r2, s_p=s_p, h=h, Lam=c_0 * (r1 - r2))


def section_run():
    t0 = time.time()
    E = 17.48
    base = setup(E)
    g_sym = base["kappa"] / (2 * np.pi)
    out = dict(E=E, g_sym=g_sym, cases={})
    worst_map, worst_H = 0.0, 0.0
    for b in (0.5, 2.0, 4.0):
        geo = geometry(E, b)
        t_unit = base["t_unit"] * geo["lam_ratio"]  # Lambda_asym / v
        st = dict(base)
        st["t_unit"] = t_unit
        st["q0"] = 2 * np.pi * 87e9 * t_unit
        st["tau"] = 0.9 * 2 * np.pi / st["q0"]
        st["eta"] = np.concatenate([2 * np.pi * F * 1e9 * t_unit, -2 * np.pi * F[[0, 5]] * 1e9 * t_unit])
        st["kappa"] = 2 * np.pi * g_sym * geo["gfac"]
        st["S"] = base["S"]  # S = eps G Lambda / q0 is unit-independent: same peak strain 1e-5
        c_h, c_0 = 1j * np.pi * np.sqrt(b), 1j * np.pi / np.sqrt(b)
        ra = integrate(st, c_h, c_0)
        rs = integrate(st, 1j * np.pi, 1j * np.pi)
        dmap = float(np.max(np.abs(ra["R"] - np.sqrt(b) * rs["R"])))
        worst_map = max(worst_map, dmap)
        # weights at 14 ns (physical) through the continuation, both runs; asym normalised with c_h (kinematic law with alpha sqrt(b))
        T = 14e-9 / t_unit
        wa, Ta = tt_weight(st, ra, T)
        ws, _ = tt_weight(st, rs, T)
        hn = h_norm(st, Ta)  # |G_T| |R0~| pi S / kappa, with the symmetric root
        H_sym = ws / hn
        H_asym = wa / (hn * b)  # |R0_asym| = sqrt(b) |R0~| and alpha_h = sqrt(b): kinematic law x b
        _, B = closed_forms(st, Ta)
        nf = F.size
        dH = float(np.max(np.abs(H_asym[:nf] / B[:nf] - 1)))
        worst_H = max(worst_H, dH)
        print(f"b = {b:g}: chi = {geo['chi_deg']:+.2f} deg, Lambda_asym/Lambda_sym = {geo['lam_ratio']:.3f}, g_asym = g ({geo['gfac']:.4f}) = {g_sym * geo['gfac']:.5f}; "
              f"mapping: max|R_asym - sqrt(b) R~| = {dmap:.1e}; W_asym / (b W~) at 14 ns: {np.max(np.abs(wa / (b * ws) - 1)):.1e}; "
              f"H_TT(b) / CF-B(g_asym) worst {dH * 100:.2f} %  [H(20)/H(100) TT {H_asym[0] / H_asym[5]:.4f}, CF {B[0] / B[5]:.4f}]  ({time.time() - t0:.0f} s)", flush=True)
        out["cases"][f"b{b:g}"] = dict(geo=geo, map_err=dmap, H_err=dH, H_asym=H_asym.tolist(), H_cf=B.tolist(), T14_units=Ta)
    out["verdict_map"] = "mapping exact (< 1e-6)" if worst_map < 1e-6 else "MAPPING WRONG"
    out["verdict_H"] = "one geometry factor" if worst_H < 0.003 else "closed form needs extension"
    print("verdict:", out["verdict_map"], "/", out["verdict_H"])
    # closed-form table over b: steady ratio and window boundaries (in ns, with Lambda_asym)
    table = {}
    for b in (0.5, 1.0, 2.0, 4.0, 10.0):
        geo = geometry(E, b)
        t_unit = base["t_unit"] * geo["lam_ratio"]
        st = dict(base)
        st["t_unit"] = t_unit
        st["q0"] = 2 * np.pi * 87e9 * t_unit
        st["tau"] = 0.9 * 2 * np.pi / st["q0"]
        st["eta"] = 2 * np.pi * F * 1e9 * t_unit
        st["kappa"] = 2 * np.pi * g_sym * geo["gfac"]
        st["S"] = base["S"]
        steady = closed_forms(st, 1e4)[1]
        r_inf = steady[0] / steady[5]
        Tg = np.geomspace(0.5, 500, 1500) * 1e-9 / t_unit
        ratio = np.array([(lambda B: B[0] / B[5])(closed_forms(st, T)[1]) for T in Tg])
        bnd = {}
        for name, target in (("T_flat", 0.9), ("T_half", 0.5 * (1 + r_inf)), ("T_90", r_inf + 0.1 * (1 - r_inf))):
            j = int(np.argmax(ratio < target))
            bnd[name] = float(Tg[j] * t_unit * 1e9)
        table[f"b{b:g}"] = dict(chi_deg=geo["chi_deg"], lam_um=base["p"]["Lambda_um"] * geo["lam_ratio"], g_asym=g_sym * geo["gfac"], H20_H100=float(r_inf), H20=float(steady[0]), H100=float(steady[5]), **bnd)
        print(f"  b = {b:4g}: chi = {geo['chi_deg']:+6.2f} deg, Lambda = {table[f'b{b:g}']['lam_um']:.2f} um, g = {g_sym * geo['gfac']:.5f}, steady H(20) = {steady[0]:.3f}, H(100) = {steady[5]:.3f}, ratio {r_inf:.3f}, "
              f"T_flat {bnd['T_flat']:.1f} / T_half {bnd['T_half']:.1f} / T_90 {bnd['T_90']:.1f} ns")
    out["table"] = table
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT.name} ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run"]:
        globals()[f"section_{sec}"]()
