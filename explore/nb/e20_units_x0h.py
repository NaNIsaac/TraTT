"""Entry 20: experiment mapping with verified X0h constants (Si, sigma polarization).

X0h (x-server.gmca.aps.anl.gov, queried 2026-10-01), symmetric Laue:
ext_L = lambda cos(theta_B) / |chi_rh| is the intensity Pendelloesung period, i.e. the
code's Lambda_P (unmodulated I = sin^2(pi L / Lambda_P)); checked by hand for Si 220 8 keV.

Key identity: alpha / alpha_ac = Lambda_ac / Lambda_P. The resolved-channel rule
(Entry 19: <= 0.05 for 3 % accuracy, <= 0.1 for 8 %) is a bound on the acoustic
wavelength alone, so the minimum drive frequency f = v_s / Lambda_ac is set by Lambda_P,
independent of the crystal thickness.
"""
import numpy as np
from scipy.special import jv

V_T = np.sqrt(79.6e9 / 2329.0)  # shear wave along [001], c44 / rho -> 5846 m/s
D = {"220": 1.9201e-10, "400": 1.3577e-10, "440": 0.9601e-10}
# (E keV, hkl): Lambda_P [um], theta_B [deg], |chi_rh|, |chi_ih|, chi_i0, mu0 [1/cm]
X0H = {
    (8, "220"): (15.203, 23.802, 0.93209e-5, 0.34424e-6, 0.35783e-6, 145.07),
    (8, "400"): (16.285, 34.802, 0.78073e-5, 0.33115e-6, 0.35783e-6, 145.07),
    (8, "440"): (15.416, 53.818, 0.59269e-5, 0.31001e-6, 0.35783e-6, 145.07),
    (17.48, "220"): (36.288, 10.644, 0.19209e-5, 0.15496e-7, 0.16059e-7, 14.226),
    (17.48, "400"): (42.659, 15.142, 0.16049e-5, 0.14952e-7, 0.16059e-7, 14.226),
    (17.48, "440"): (54.336, 21.679, 0.12130e-5, 0.13921e-7, 0.16059e-7, 14.226),
    (30, "220"): (63.339, 6.1781, 0.64870e-6, 0.17580e-8, 0.18195e-8, 2.7663),
    (30, "440"): (98.759, 12.430, 0.40866e-6, 0.15857e-8, 0.18195e-8, 2.7663),
}
S = 1.5
print(f"v_T = {V_T:.0f} m/s, S = {S}, L = 3 Lambda_P (alpha = 3, alpha_1 = {3*jv(1,S):.2f})")
print("E keV hkl | Lambda_P um | L um | t_tr ns | tau_1 ns | f_min GHz (0.05 / 0.1) | strain at 0.1 | "
      "kappa | a''/a' | contrast t=1 | exp(-mu L/cos)")
for (E, hkl), (lp, th, xr, xi, xi0, mu) in X0H.items():
    Lp = lp * 1e-6
    L = 3 * Lp
    t_tr = L / V_T
    tau1 = Lp / (V_T * jv(1, S))
    f05, f10 = V_T / (0.05 * Lp), V_T / (0.1 * Lp)
    strain = S * D[hkl] / (0.1 * Lp)          # K u0 = S d / Lambda_ac
    a_im = 3 * xi / xr                         # alpha'' with alpha' = 3
    contrast = 1 / np.cosh(2 * np.pi * a_im * jv(1, S))
    trans = np.exp(-mu * 100 * L / np.cos(np.radians(th)))
    print(f"{E:5g} {hkl} | {lp:7.2f} | {L*1e6:6.1f} | {t_tr*1e9:6.2f} | {tau1*1e9:6.2f} | "
          f"{f05/1e9:5.2f} / {f10/1e9:5.2f} | {strain:.1e} | {xi/xi0:.3f} | {xi/xr:.4f} | "
          f"{contrast:.3f} | {trans:.3f}")
