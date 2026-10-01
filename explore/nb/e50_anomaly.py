"""Entry 50: anomaly ranking over the phase-averaged phase diagram (cache_e43_ph4).

Two reference models for the order-window populations p_n:
  IND   independent channels, coupling pi alpha |J_n(S)|
  G3    independent channels with the third-order effective coupling
        |g_n| = pi alpha |J_{-n}(S) - (R/2)^2 C_n(S)|   (Entry 45), level shifts ignored
Scores per point (alpha_ac >= 2, so that windows are channels at small R):
  E_ind, E_g3       population errors vs the two models
  gain = E_ind - E_g3 (what the indirect coupling explains)
  dark  = min_n p_n / p_n^G3 over orders with |J_n| > 0.2   (suppressed despite large coupling)
  bright= max_n p_n / max(p) over orders with |J_n| < 0.05  (bright near a zero)
Prints the top cases by E_g3 (unexplained even with indirect coupling), by darkness and brightness,
and the regime split of E_g3.
"""
import glob
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e45_zero import C_n  # noqa: E402


def window_laue(a, q):
    """Integral of the thickness-1 Laue line of coupling a over |d| < q/2."""
    d = np.linspace(-q / 2, q / 2, 2001)
    s2 = a * a + d * d / 4
    return np.trapezoid(a * a * np.sin(np.sqrt(s2)) ** 2 / np.where(s2 > 0, s2, 1), d) if a > 0 else 0.0


if __name__ == "__main__":
    rows = [json.loads(open(f).read()) for f in glob.glob(str(HERE / "cache_e43_ph4" / "*.json"))]
    Ccache = {}
    out = []
    for r in rows:
        S, alpha, aac = r["S"], r["alpha"], r["aac"]
        if S == 0 or aac < 2:
            continue
        q, R = 2 * np.pi * aac, alpha / aac
        ns, p, pi = np.array(r["ns"]), np.array(r["p"]), np.array(r["p_ind"])
        g3 = []
        for n in ns:
            key = (int(n), S)
            if key not in Ccache:
                Ccache[key] = C_n(int(n), S, K=25)
            g3.append(window_laue(np.pi * alpha * abs(jv(-n, S) - (R / 2) ** 2 * Ccache[key]), q))
        g3 = np.array(g3)
        E_ind = np.abs(p - pi).sum() / pi.sum()
        E_g3 = np.abs(p - g3).sum() / g3.sum()
        big = np.abs(jv(ns, S)) > 0.2
        small = (np.abs(jv(ns, S)) < 0.05) & (np.abs(ns) <= S + 2)
        dark = float(np.min(p[big] / g3[big])) if big.any() else np.nan
        nd = int(ns[big][np.argmin(p[big] / g3[big])]) if big.any() else 0
        bright = float(np.max(p[small] / p.max())) if small.any() else np.nan
        out.append(dict(alpha=alpha, aac=aac, S=S, R=R, E_ind=E_ind, E_g3=E_g3, dark=dark, n_dark=nd, bright=bright))
    json.dump(out, open(HERE / "e50_anomaly.json", "w"))
    R = np.array([o["R"] for o in out])
    Ei, Eg = np.array([o["E_ind"] for o in out]), np.array([o["E_g3"] for o in out])
    print("median E_ind / E_g3 by R:")
    for lo, hi in [(0, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1), (1, 2), (2, 30)]:
        m = (R > lo) & (R <= hi)
        print(f"  R ({lo},{hi}]: {np.median(Ei[m]):.4f} / {np.median(Eg[m]):.4f}  (n={m.sum()})")
    for key, rev, title in (("E_g3", True, "unexplained even with third-order coupling"),
                            ("dark", False, "darkest strong orders (vs G3)"), ("bright", True, "brightest near-zero orders")):
        print(f"== top 10 {title} ==")
        good = [o for o in out if np.isfinite(o[key])]
        for o in sorted(good, key=lambda o: o[key], reverse=rev)[:10]:
            print(f"  alpha={o['alpha']:5.2f} aac={o['aac']:4.1f} S={o['S']:3.1f} R={o['R']:5.2f} | E_ind {o['E_ind']:.3f} "
                  f"E_g3 {o['E_g3']:.3f} | dark {o['dark']:.3f} (n={o['n_dark']}) | bright {o['bright']:.3f}")
