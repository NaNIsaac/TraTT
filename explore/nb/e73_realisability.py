"""Entry 73 (T1c): what S = 0.5-6 means physically for Si 440 (and 400, 220) at 1-2 GHz.

Pure constants, no TT. Shear wave along the surface normal, polarised along h (symmetric Laue, Entry 6/20),
v_T = sqrt(c44 / rho) = 5846 m/s (textbook c44 = 79.6 GPa, rho = 2329 kg/m^3; not re-verified).
    u0        = S d_hkl / (2 pi)                       absolute displacement amplitude
    eps       = K u0 = S d_hkl / Lambda_ac             shear strain amplitude, Lambda_ac = v_T / f
    sigma     = c44 eps                                 shear stress amplitude
    w         = (1/2) c44 eps^2                         time-averaged acoustic energy density (travelling wave)
    I         = v_T w                                   acoustic intensity
    eps2/eps  = (beta/2) k eps L                        second-harmonic strain ratio after L (beta ~ 5, estimate)
    R         = Lambda_ac / Lambda_P                    independent-satellite validity (<= 0.1, Entry 19/54)
    beta_qs   = 2 pi f L / c                            quasi-static X-ray condition (Entry 62)
L = 3 Lambda_P as in Entry 20. The "large displacement vs large strain" distinction is u0 / d against eps.
Usage: python e73_realisability.py
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e20_units_x0h import V_T, X0H, D  # noqa: E402

C44 = 79.6e9
BETA_NL = 5.0  # acoustic nonlinearity parameter, order of magnitude (not verified)
C = 2.998e8
S_GRID = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0]
CASES = [(17.48, "440"), (30, "440"), (17.48, "400"), (17.48, "220")]
F_GHZ = [1.0, 2.0]

rows = []
print(
    "E keV hkl  f GHz | Lambda_ac um  R=Lac/LP  beta_qs | S | u0 pm  u0/d | strain | stress MPa | w J/m3 | I W/cm2 | eps2/eps over L"
)
for E, hkl in CASES:
    lp = X0H[(E, hkl)][0] * 1e-6
    L = 3 * lp
    d = D[hkl]
    for f in F_GHZ:
        lac = V_T / (f * 1e9)
        R = lac / lp
        bqs = 2 * np.pi * f * 1e9 * L / C
        k = 2 * np.pi / lac
        for S in S_GRID:
            u0 = S * d / (2 * np.pi)
            eps = S * d / lac
            sig = C44 * eps
            w = 0.5 * C44 * eps**2
            inten = V_T * w
            h2 = 0.5 * BETA_NL * k * eps * L
            rows.append(
                dict(
                    E=E,
                    hkl=hkl,
                    f=f,
                    lac_um=lac * 1e6,
                    R=R,
                    beta_qs=bqs,
                    S=S,
                    u0_pm=u0 * 1e12,
                    u0_over_d=u0 / d,
                    strain=eps,
                    stress_MPa=sig / 1e6,
                    w=w,
                    I_W_cm2=inten / 1e4,
                    h2=h2,
                )
            )
            print(
                f"{E:5g} {hkl} {f:3g}   | {lac * 1e6:6.2f}  {R:6.3f}  {bqs:.1e} | {S:3g} | {u0 * 1e12:5.1f} {u0 / d:6.3f} | "
                f"{eps:.1e} | {sig / 1e6:6.2f} | {w:7.3f} | {inten / 1e4:7.4f} | {h2:.3f}"
            )
json.dump(rows, open(HERE / "e73_realisability.json", "w"), indent=1)
# transducer rise: a resonator of quality factor Q rings up over ~Q / f
print("\nrise time Q/f vs t_tr = L / v_T:")
for E, hkl in CASES[:2]:
    lp = X0H[(E, hkl)][0] * 1e-6
    t_tr = 3 * lp / V_T
    for f in F_GHZ:
        print(
            f"  {E:5g} {hkl} {f:g} GHz: t_tr = {t_tr * 1e9:.1f} ns; Q = 10 / 30 / 100 -> rise/t_tr = "
            + " / ".join(f"{Q / (f * 1e9) / t_tr:.2f}" for Q in (10, 30, 100))
        )
