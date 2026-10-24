"""Entry 85 (T2o): one GHz-train case, Si 440 at 17.48 keV, 1-2 GHz, in one table. No new physics: the formulas and
constants of Entries 20 (X0h), 73 (elastic numbers), 75 (window law, CF-B), 82-84 (beam conditions).

Geometry of the poster: symmetric Laue, L = 3 Lambda_P (alpha = 3), bulk shear wave along the surface normal,
v_T = 5846 m/s; the Bragg numbers (window law) are for the same reflection in symmetric Bragg.
Usage: python e85_si440_case.py
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e75_window_cf import closed_forms  # noqa: E402

OUT = HERE / "e85_si440_case.json"
# X0h, Si 440, 17.48 keV, sigma (e20_units_x0h.py): Lambda_P_Laue um, theta_B deg, |chi_rh|, |chi_ih|, chi_i0, mu0 /cm
LAM_P, TH_B, XRH, XIH, XI0, MU0 = 54.336, 21.679, 0.12130e-5, 0.13921e-7, 0.16059e-7, 14.226
D440, HC_E = 0.9601e-10, 12.398419843 / 17.48 * 1e-10  # d spacing, wavelength (m)
V_T, C44, RHO = 5846.0, 79.6e9, 2329.0
# Bragg counterpart: in symmetric Bragg 440 h is along the surface normal, so a shear wave along the normal has
# h.u = 0 (no coupling; physics cross-review of Entry 85). The Bragg window numbers use a longitudinal wave
# along [110], v_L = sqrt((c11 + c12 + 2 c44) / 2 rho) = 9136 m/s (c11 = 165.7, c12 = 63.9 GPa).
V_L110 = 9136.0
ALPHA = 3.0
L = ALPHA * LAM_P * 1e-6
F_GHZ = (1.0, 1.5, 2.0)
S_LIST = (0.01, 0.02, 0.05, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0)


def main():
    th = np.radians(TH_B)
    g = XI0 / XRH
    darwin_half = XRH / np.sin(2 * th)  # rad
    meV_per_hw = darwin_half / np.tan(th) * 17.48e6
    lam_bragg = HC_E * np.sin(th) / XRH  # symmetric Bragg extinction length (m)
    t_tr = L / V_T
    out = dict(g=g, darwin_half_urad=darwin_half * 1e6, meV_per_hw=meV_per_hw, lam_bragg_um=lam_bragg * 1e6, t_tr_ns=t_tr * 1e9, rows=[], window={})
    print(f"Si 440, 17.48 keV: Lambda_P(Laue) = {LAM_P} um, L = 3 Lambda_P = {L * 1e6:.0f} um, t_tr = {t_tr * 1e9:.1f} ns, g = chi_i0/|chi_rh| = {g:.4f}, "
          f"Darwin half-width {darwin_half * 1e6:.2f} urad = {meV_per_hw:.0f} meV, Lambda(Bragg) = {lam_bragg * 1e6:.1f} um, f_D(Bragg, longitudinal v_L[110] = {V_L110:.0f} m/s) = {V_L110 / lam_bragg / 1e9:.3f} GHz")
    print("\n f GHz |  S   | u0 pm | strain  | I W/cm2 | Lambda_ac um | alpha_ac |   R   | nu=pi a |J1(S)| | tau_1/t_tr | tau_1 ns | beta_qs | beam(P2 centre) 0.3 w_1 urad | beam(ring, Bragg longitudinal) g+1/2x urad")
    for f in F_GHZ:
        lam_ac = V_T / (f * 1e9)
        aac = L / lam_ac
        R = lam_ac / (LAM_P * 1e-6)
        bqs = 2 * np.pi * f * 1e9 * L / 2.998e8
        x_bragg = f * 1e9 * lam_bragg / V_L110  # Bragg counterpart, longitudinal wave
        ring_hw = g + 1 / (2 * x_bragg)
        for S in S_LIST:
            u0 = S * D440 / (2 * np.pi)
            eps = S * D440 / lam_ac
            I = V_T * 0.5 * C44 * eps**2 / 1e4
            an = ALPHA * abs(jv(1, S))
            nu = np.pi * an
            tau1 = 0.5 / an
            w1_hw = abs(jv(1, S))
            row = dict(f_GHz=f, S=S, u0_pm=u0 * 1e12, strain=eps, I_Wcm2=I, lam_ac_um=lam_ac * 1e6, alpha_ac=aac, R=R, nu=nu, tau1_ttr=tau1, tau1_ns=tau1 * t_tr * 1e9,
                       beta_qs=bqs, beam_centre_urad=0.3 * w1_hw * darwin_half * 1e6, beam_ring_urad=ring_hw * darwin_half * 1e6, x_bragg=x_bragg)
            out["rows"].append(row)
            print(f" {f:5g} | {S:4g} | {u0 * 1e12:5.1f} | {eps:.1e} | {I:7.2f} | {lam_ac * 1e6:12.2f} | {aac:8.1f} | {R:.3f} | {nu:14.3f} | {tau1:10.3f} | {tau1 * t_tr * 1e9:8.1f} | {bqs:.1e} | {0.3 * w1_hw * darwin_half * 1e6:28.2f} | {ring_hw * darwin_half * 1e6:.2f}")
    # window law for the same reflection in symmetric Bragg (CF-B, Entry 75), ratio H(1 GHz)/H(2 GHz)
    t_unit = lam_bragg / V_L110
    st = dict(t_unit=t_unit, kappa=2 * np.pi * g, S=1e-5 * 2 * np.pi / D440 * lam_bragg / (2 * np.pi * 87e9 * t_unit) / 1000, q0=2 * np.pi * 87e9 * t_unit, tau=0.9 * 2 * np.pi / (2 * np.pi * 87e9 * t_unit),
              eta=2 * np.pi * np.array([1e9, 2e9]) * t_unit)
    steady = closed_forms(st, 1e4)[1]
    r_inf = steady[0] / steady[1]
    Tg = np.geomspace(1, 2000, 2000) * 1e-9 / t_unit
    ratio = np.array([(lambda B: B[0] / B[1])(closed_forms(st, T)[1]) for T in Tg])
    for name, target in (("T_flat", 0.9), ("T_half", 0.5 * (1 + r_inf)), ("T_90", r_inf + 0.1 * (1 - r_inf))):
        j = int(np.argmax(ratio < target))
        out["window"][name] = float(Tg[j] * t_unit * 1e9)
    out["window"]["H1_over_H2_steady"] = float(r_inf)
    out["window"]["H_steady_1GHz"] = float(steady[0])
    out["window"]["H_steady_2GHz"] = float(steady[1])
    out["window"]["kappa_depth_um"] = float(lam_bragg * 1e6 / (2 * np.pi * g))
    print(f"\nwindow law, symmetric Bragg Si 440 / 17.48 keV, longitudinal [110] wave (a shear wave has h.u = 0 there): x = {st['eta'][0] / (2 * np.pi):.2f} (1 GHz), {st['eta'][1] / (2 * np.pi):.2f} (2 GHz); steady H = {steady[0]:.3f} / {steady[1]:.3f}, ratio {r_inf:.3f}; "
          f"T_flat {out['window']['T_flat']:.0f} ns, T_half {out['window']['T_half']:.0f} ns, T_90 {out['window']['T_90']:.0f} ns; 1/kappa = {out['window']['kappa_depth_um']:.0f} um = {out['window']['kappa_depth_um'] * 1e-6 / V_L110 * 1e9:.0f} ns")
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT.name}")


if __name__ == "__main__":
    main()
