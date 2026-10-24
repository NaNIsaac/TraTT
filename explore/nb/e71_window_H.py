"""Entry 71: does the Entry 65 crystal factor H(f) survive a finite measurement window?

Entry 65 computed H for a steady single mode filling the crystal, where the satellite integrates over the
absorption depth 1/kappa (16 Lambda = 185 um at Si 400, 17.48 keV). Entry 70 found the 8 keV and 17.48 keV
packet responses identical over 140 ps and 1400 ps windows: in that time the packet is only 0.1-1 Lambda
deep and absorption has not acted. So the factor-2 bias quoted for the Chalise et al. band may not apply
to a 140 ps measurement. Test: the measured transfer function
    H_meas(f) = W(f) / (|r0(f)| |phi^(f)|),  phi^ = eps^ / omega (Entry 59 sec. 1), normalised at 100 GHz,
from the full semi-infinite TT run of the Entry 66 packet in the linear regime (S/100), for windows
140 ps, 1.4 ns, 14 ns, against the steady closed form H(f) of Entry 65 at the same normalisation.
Pre-registered rule: Entry 65's number for their case stands only if H_meas(20)/H_meas(100) at 140 ps is
within 20 % of the steady value 0.51 (17.48 keV). If it is near 1, Entry 65 section 2 is wrong for
short windows and must be corrected.
Usage: python e71_window_H.py
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e65_si400 import H, params, root  # noqa: E402
from e66_harmonic import surface_R  # noqa: E402

V_L = 8433.0
F = np.array([20.0, 30.0, 40.0, 60.0, 87.0, 100.0, 150.0])
OUT = HERE / "e71_window_H.json"


def run(E, T_ps, ncyc=0.9, f0=87.0):
    p = params(E)
    t_unit = p["Lambda_um"] * 1e-6 / V_L
    g_lambda = 2 * np.pi / 1.3577e-10 * p["Lambda_um"] * 1e-6
    q0 = 2 * np.pi * f0 * 1e9 * t_unit
    eta = 2 * np.pi * F * 1e9 * t_unit
    kappa = 2 * np.pi * p["g"]
    S = 1e-3 * g_lambda / q0 / 100.0  # linear regime
    tau_d = ncyc * 2 * np.pi / q0
    T = T_ps * 1e-12 / t_unit
    g = lambda s: np.where(s > 0, np.exp(-np.clip(s, 0, None) / tau_d) * np.sin(q0 * s), 0.0)  # noqa: E731
    per_unit = 64 * eta.max() / (2 * np.pi)
    t, R = surface_R(1.0, kappa, eta, lambda s: S * g(s), T, per_unit)
    inten = np.abs(R) ** 2
    inten -= inten.mean(axis=0)
    w = np.full(t.size, t[1] - t[0])
    w[[0, -1]] /= 2
    W = np.abs(np.sum(w[:, None] * inten * np.exp(1j * eta[None, :] * t[:, None]), axis=0))
    s = np.linspace(0, T, 200001)
    eps = np.gradient(S * g(s), s)
    E_spec = np.abs([np.trapezoid(eps * np.exp(1j * o * s), s) for o in eta])
    r0 = np.abs([root(1.0, e, kappa) for e in eta])
    Hm = W * eta / (r0 * E_spec)  # W / (|r0| |phi^|), phi^ = eps^ / omega (Entry 59 sec. 1)
    i100 = list(F).index(100.0)
    return Hm / Hm[i100], p


if __name__ == "__main__":
    out = {"f_GHz": F.tolist()}
    for E in (17.48, 8.0):
        p = params(E)
        Hs = H(F / p["f_D_GHz"], p["g"])
        Hs = Hs / Hs[list(F).index(100.0)]
        print(
            f"Si 400, {E:g} keV (g = {p['g']:.4f}): steady H(f)/H(100) = "
            + " ".join(f"{f:g}:{h:.3f}" for f, h in zip(F, Hs, strict=True))
        )
        out[f"E{E}_steady"] = Hs.tolist()
        for T_ps in (140.0, 1400.0, 14000.0):
            Hm, _ = run(E, T_ps)
            print(
                f"   window {T_ps:6g} ps: H_meas/H_meas(100) = "
                + " ".join(f"{f:g}:{h:.3f}" for f, h in zip(F, Hm, strict=True)),
                flush=True,
            )
            out[f"E{E}_T{T_ps:g}"] = Hm.tolist()
    json.dump(out, open(OUT, "w"), indent=1)
