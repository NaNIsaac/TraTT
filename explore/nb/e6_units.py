"""Entry 6: order-of-magnitude mapping to a real Laue experiment.

Inputs are textbook-level numbers that MUST be re-checked against a dynamical
diffraction calculator (e.g. XOP / Stepanov X0h) before they go into the paper.
Convention of this code: unmodulated Laue I = sin^2(pi alpha) -> alpha = L / Lambda_P,
Lambda_P the intensity Pendelloesung period.
"""
import numpy as np
from scipy.special import jv

LAMBDA_P = 15e-6   # m, Si(220) symmetric Laue near 8 keV, sigma pol. -- TO VERIFY
V_T = 5.84e3       # m/s, Si transverse [100] -- shear wave needed so h.u != 0 (handoff 7.4)
D220 = 1.920e-10   # m

for L in (50e-6, 100e-6, 200e-6):
    alpha = L / LAMBDA_P
    t_tr = L / V_T
    aac = max(20.0, 5 * alpha)           # keep channels resolved: alpha_ac >= 5 alpha
    lam_ac = L / aac
    f = V_T / lam_ac
    for S in (0.5, 1.5):
        an = alpha * jv(1, S)
        u0 = S * D220 / (2 * np.pi)
        strain = 2 * np.pi / lam_ac * u0
        print(f"L={L*1e6:5.0f}um alpha={alpha:5.2f} t_tr={t_tr*1e9:5.1f}ns alpha_ac={aac:4.0f} "
              f"f={f/1e9:5.2f}GHz | S={S}: u0={u0*1e12:5.1f}pm strain={strain:.1e} "
              f"alpha_1={an:4.2f} period t_tr/alpha_1={t_tr/an*1e9:5.2f}ns "
              f"first max {t_tr/(2*an)*1e9:5.2f}ns")
