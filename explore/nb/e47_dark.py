"""Entry 47: a dark main line with a strong direct coupling (alpha 2, alpha_ac 2, S 1, n = 0).

J_0(1) = 0.765 and pi alpha J_0 = 4.8, yet the eta-window |eta| < q/2 holds ~1 % of the
independent-channel value at every acoustic phase (Entry 43 cache, 8-phase check).
Diagnostics (phase-averaged over 4 phases unless stated):
  (1) sensitivity to S, alpha_ac, alpha;
  (2) path removal in the Floquet basis;
  (3) the exact two-band gap at eta = 0 (monodromy, Entry 46);
  (4) where the weight is: spectrum near eta = 0 vs the independent line.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.special import jv  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e43_phase import laue_line  # noqa: E402
from e46_bands import bands  # noqa: E402
from tratt.core import required_nz, solve_laue  # noqa: E402
from tratt.drive import Drive  # noqa: E402
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402


def main_window(alpha, aac, S, nph=4, mask=None, deta=0.1):
    q = 2 * np.pi * aac
    eta = np.arange(-q / 2, q / 2, deta)
    if mask is None:
        nz = required_nz(alpha, aac, S=S, eta_max=q)
        I = np.mean([np.abs(solve_laue(alpha, eta, phi=Drive(S=S, alpha_ac=aac).steady_phi(k / (nph * aac)),
                                       n_z=nz)[1]) ** 2 for k in range(nph)], axis=0)
    else:
        K = auto_K(S) + 4 + int(np.ceil(4 * np.pi * alpha / q))
        I = np.mean([FloquetLaue(alpha, aac, S, K, t=k / (nph * aac), mask=mask).rocking(eta)
                     for k in range(nph)], axis=0)
    ind = laue_line(eta, np.pi * alpha * abs(jv(0, S)))
    return np.trapezoid(I, eta), np.trapezoid(ind, eta), eta, I, ind


if __name__ == "__main__":
    a0, q0, S0 = 2.0, 2.0, 1.0
    p, pi, eta, I, ind = main_window(a0, q0, S0)
    print(f"base: window integral {p:.4f} vs independent {pi:.3f} (ratio {p/pi:.4f}); max I {I.max():.3f} vs {ind.max():.3f}")
    print("(1) sensitivity, ratio p/p_ind:")
    print("   S:     " + "  ".join(f"{S:.2f}:{(lambda r: r[0]/r[1])(main_window(a0, q0, S)):.3f}" for S in np.linspace(0.5, 1.6, 12)))
    print("   aac:   " + "  ".join(f"{q:.2f}:{(lambda r: r[0]/r[1])(main_window(a0, q, S0)):.3f}" for q in (1.5, 1.75, 1.9, 2.0, 2.1, 2.25, 2.5, 3.0)))
    print("   alpha: " + "  ".join(f"{a:.2f}:{(lambda r: r[0]/r[1])(main_window(a, q0, S0)):.3f}" for a in (1.0, 1.5, 1.8, 2.0, 2.2, 2.5, 3.0)))
    print("(2) path removal (Floquet, phase-avg):")
    for lab, mk in (("full", lambda k, m: True), ("no forward sidebands (STAR)", lambda k, m: k == 0),
                    ("|k| <= 1", lambda k, m: abs(k) <= 1), ("no b_{+-1}", lambda k, m: abs(m) != 1),
                    ("no b_{+-2}", lambda k, m: abs(m) != 2), ("only a_0, b_0 (two-beam)", lambda k, m: k == 0 and m == 0),
                    ("only a_0, b_{0,+-1}", lambda k, m: k == 0 and abs(m) <= 1)):
        r = main_window(a0, q0, S0, mask=mk)
        print(f"   {lab:28s}: p/p_ind {r[0]/r[1]:.4f}")
    q = 2 * np.pi * q0
    e = np.array([0.0, 0.25 * q, 0.5 * q])
    lam, gap, mix, _ = bands(a0, q0, S0, e)
    print(f"(3) two-band gap at eta = 0, q/4, q/2: {np.round(gap, 3)} (independent 2 pi alpha J0 = {2*np.pi*a0*jv(0,S0):.2f}, q/2 = {q/2:.2f});"
          f" diffracted weight of the two Bloch modes at eta=0: {np.round(mix[0], 3)}")
    pw, _, ew, Iw, _ = main_window(a0, q0, S0, deta=0.1)
    big = (lambda e: np.arange(-1.5 * q, 1.5 * q, 0.1))(0)
    nz = required_nz(a0, q0, S=S0, eta_max=1.5 * q)
    Ib = np.mean([np.abs(solve_laue(a0, big, phi=Drive(S=S0, alpha_ac=q0).steady_phi(k / (4 * q0)), n_z=nz)[1]) ** 2
                  for k in range(4)], axis=0)
    print(f"(4) weight in |eta|<q/2: {np.trapezoid(np.where(np.abs(big)<q/2, Ib, 0), big):.3f}; in q/2<|eta|<3q/2: "
          f"{np.trapezoid(np.where((np.abs(big)>=q/2), Ib, 0), big):.3f}; independent total over the same range "
          f"{np.trapezoid(sum(laue_line(big - n*q, np.pi*a0*abs(jv(n,S0))) for n in range(-4,5)), big):.3f}")
    np.savez(HERE / "e47_dark.npz", eta=big, I=Ib, q=q)
