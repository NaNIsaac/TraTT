"""Entry 53: is the 4-phase acoustic average converged, and does the third-order residual coupling
at a Bessel zero dominate the physical exit intensity?

Definitions (existing code):
  Floquet H[a_k,b_m] = pi alpha J_{k+m}(S) e^{-i(k+m) q t}   (src/tratt/floquet.py:63), x(0) = a_0,
  Dh(zeta) = sum_m b_m(zeta) e^{-i m q zeta}                 (floquet.py:94)
  TT steady phase phi(zeta) = S sin(q zeta - q t)             (src/tratt/drive.py, steady_phi)
  acoustic phase phi_ac = q t; Entry 43 averages intensity over t = k / (N alpha_ac), i.e.
  phi_ac = 2 pi k / N, k = 0..N-1 (explore/nb/e43_phase.py:62-64), N = 4.
Result derived in notebook Entry 53: with D = diag(e^{i k phi_ac} on a_k, e^{-i m phi_ac} on b_m),
D H(phi_ac) D^dagger = H(0) and D a_0 = a_0, so b_m(1; phi_ac) = e^{i m phi_ac} b_m(1; 0) and
    A_h(phi_ac) = Dh(1) = sum_m c_m e^{i m phi_ac},   c_m = b_m(1; 0) e^{-i m q},
    <|A_h|^2>_{phi_ac uniform} = sum_m |b_m(1; 0)|^2,
while an N-point uniform average keeps the cross terms with m - m' = 0 (mod N).
Usage: python e53_phaseavg.py [check] [eps] [map]
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.special import jn_zeros, jv  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e52_heff import B_n  # noqa: E402
from tratt.core import required_nz, solve_laue  # noqa: E402
from tratt.drive import Drive  # noqa: E402
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402

J11 = jn_zeros(1, 1)[0]


def components(alpha, aac, S, eta, K=None):
    """b_m(1; phi_ac = 0) for every eta: shape (n_eta, 2K+1), and the m labels."""
    q = 2 * np.pi * aac
    K = K or auto_K(S) + 4 + int(np.ceil(4 * np.pi * alpha / q)) + int(np.ceil(np.abs(np.asarray(eta)).max() / q))
    fl = FloquetLaue(alpha, aac, S, K)
    x = fl.propagate(np.atleast_1d(eta))
    return x[:, fl.n:], fl.ms, q


def phase_avg_N(b, ms, q, N):
    """N-point uniform average of |sum_m c_m e^{i m phi}|^2 with c_m = b_m e^{-i m q}."""
    c = b * np.exp(-1j * ms * q)[None, :]
    ph = 2 * np.pi * np.arange(N) / N
    A = c @ np.exp(1j * np.outer(ms, ph))
    return np.mean(np.abs(A) ** 2, axis=1)


def section_check():
    print("== reproduction and convergence: alpha = 1, S = j11, channel n = 1 (target b_-1) ==")
    alpha, n = 1.0, 1
    rows = []
    for aac in (20.0, 20.25, 20.5):
        q = 2 * np.pi * aac
        R = alpha / aac
        dstar = -np.pi * alpha * R * B_n(n, J11)
        for lab, e0 in (("nominal eta = q", n * q), ("corrected q + delta*", n * q + dstar)):
            eta = np.array([e0])
            b, ms, _ = components(alpha, aac, J11, eta)
            b2, ms2, _ = components(alpha, aac, J11, eta, K=auto_K(J11) + 14)
            tgt = abs(b[0, list(ms).index(-n)]) ** 2
            exact = float(np.sum(np.abs(b[0]) ** 2))
            exact2 = float(np.sum(np.abs(b2[0]) ** 2))
            avgs = {N: float(phase_avg_N(b, ms, q, N)[0]) for N in (4, 8, 16, 32, 64)}
            # direct TT at representative phases, compared with the Floquet amplitude
            tt = []
            for k in (0, 1, 3):
                ph = 2 * np.pi * k / 4
                nz = 2 * required_nz(alpha, aac, S=J11, eta_max=e0)
                Dh = solve_laue(alpha, eta, phi=Drive(S=J11, alpha_ac=aac).steady_phi(ph / q), n_z=nz)[1]
                Af = float(phase_avg_N(b, ms, q, 1)[0]) if k == 0 else None
                c = b[0] * np.exp(-1j * ms * q)
                tt.append((abs(Dh[0]) ** 2, abs(np.sum(c * np.exp(1j * ms * ph))) ** 2))
            rows.append(dict(aac=aac, point=lab, target=tgt, exact=exact, exact_Kplus=exact2, avgs=avgs, tt=tt))
            print(f" aac={aac:6.2f} {lab:22s}: target |b_-1|^2 {tgt:.4e} | exact avg sum|b_m|^2 {exact:.4e} "
                  f"(K+ {exact2:.4e}) | N-phase: " + " ".join(f"{N}:{v:.4e}" for N, v in avgs.items()) +
                  " | TT vs Floquet at phi=0,pi/2,3pi/2: " + " ".join(f"{a:.4e}/{f:.4e}" for a, f in tt), flush=True)
    json.dump(rows, open(HERE / "e53_check.json", "w"), indent=1)


def section_eps():
    print("== eps scan: alpha = 1, S = j11, aac = 20 + eps ==")
    alpha, n = 1.0, 1
    eps = np.unique(np.concatenate([np.linspace(-0.5, 0.5, 41), np.linspace(-0.02, 0.02, 41),
                                    np.linspace(-0.002, 0.002, 21)]))
    out = []
    for e in eps:
        aac = 20.0 + e
        q = 2 * np.pi * aac
        R = alpha / aac
        dstar = -np.pi * alpha * R * B_n(n, J11)
        win = np.linspace(-0.45 * q, 0.45 * q, 1801)    # Entry 45 window, |eta - nq| < 0.45 q
        etas = np.concatenate([[n * q, n * q + dstar], n * q + win])
        b, ms, _ = components(alpha, aac, J11, etas)
        P = np.abs(b) ** 2
        it = list(ms).index(-n)
        tgt, tot = P[:, it], P.sum(1)
        w = win
        out.append(dict(eps=float(e), tgt_nom=float(tgt[0]), tot_nom=float(tot[0]), tgt_cor=float(tgt[1]),
                        tot_cor=float(tot[1]), tgt_int=float(np.trapezoid(tgt[2:], w)), tot_int=float(np.trapezoid(tot[2:], w)),
                        tgt_max=float(tgt[2:].max()), tot_at_tgtmax=float(tot[2:][np.argmax(tgt[2:])]),
                        dstar=float(dstar)))
    json.dump(out, open(HERE / "e53_eps.json", "w"), indent=1)
    for o in out[:: max(1, len(out) // 25)]:
        print(f" eps={o['eps']:+.4f}: nominal tgt {o['tgt_nom']:.3e} tot {o['tot_nom']:.3e} frac {o['tgt_nom']/o['tot_nom']:.3f} |"
              f" corrected tgt {o['tgt_cor']:.3e} tot {o['tot_cor']:.3e} frac {o['tgt_cor']/o['tot_cor']:.3f} |"
              f" window integral frac {o['tgt_int']/o['tot_int']:.4f}", flush=True)
    for key, lab in (("nom", "nominal"), ("cor", "corrected")):
        f = np.array([o[f"tgt_{key}"] / o[f"tot_{key}"] for o in out])
        e = np.array([o["eps"] for o in out])
        for thr in (0.5, 0.9):
            m = f >= thr
            print(f"  {lab}: target >= {thr:.0%} of exit intensity for eps in "
                  + (f"[{e[m].min():+.4f}, {e[m].max():+.4f}] ({m.sum()} of {m.size} points; contiguous: "
                     f"{bool(np.all(np.diff(np.flatnonzero(m)) == 1))})" if m.any() else "none"))


if __name__ == "__main__":
    for s in (sys.argv[1:] or ["check", "eps"]):
        globals()[f"section_{s}"]()
