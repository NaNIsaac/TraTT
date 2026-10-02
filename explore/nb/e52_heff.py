"""Entry 52: effective 2x2 Hamiltonian of the resonant subspace {a_0, b_-n}, derived and tested.

Exact Floquet Hamiltonian (gauge t = 0): H[a_k,a_k] = -k q, H[b_m,b_m] = eta + m q,
H[a_k,b_m] = V J_{k+m}(S), V = pi alpha, V / q = R / 2.  P = {a_0, b_-n}, eta = n q + delta.

Quasi-degenerate (Schrieffer-Wolff) expansion, derived in notebook Entry 52:
  diagonal, 2nd order:   Delta_a = -V (R/2) B_n,   Delta_b = +V (R/2) B_n,
                         B_n = sum_{k != 0} J_{k-n}^2 / k
      -> resonance at delta* = Delta_a - Delta_b = -pi alpha R B_n   (B_0 = 0)
  off-diagonal, 2nd order: identically 0 (bipartite H)
  off-diagonal, 3rd order: g = V [ J_{-n} (1 - W_n (R/2)^2 / 2) - (R/2)^2 C_n ]
      C_n = sum_{k != 0, m != -n} J_{k-n} J_{k+m} J_m / (k (m + n))
      W_n = sum_{m != -n} J_m^2 / (m + n)^2 + sum_{k != 0} J_{k-n}^2 / k^2
  diagonal, 3rd order: identically 0.
Exact reference: eigenvalues of the truncated Floquet matrix; near eta = n q the two branches
closest to 0 have minimum splitting 2 |g| at eta = n q + delta* (all orders included).
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402
from scipy.special import jn_zeros, jv, jvp  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402

KS = 40


def B_n(n, S):
    return sum(jv(k - n, S) ** 2 / k for k in range(-KS, KS + 1) if k != 0)


def C_n(n, S):
    return sum(jv(k - n, S) * jv(k + m, S) * jv(m, S) / (k * (m + n))
               for k in range(-KS, KS + 1) if k != 0 for m in range(-KS, KS + 1) if m != -n)


def W_n(n, S):
    return (sum(jv(m, S) ** 2 / (m + n) ** 2 for m in range(-KS, KS + 1) if m != -n)
            + sum(jv(k - n, S) ** 2 / k ** 2 for k in range(-KS, KS + 1) if k != 0))


def predictions(alpha, aac, S, n):
    V, r = np.pi * alpha, alpha / (2 * aac)
    g1 = V * jv(-n, S)
    g3C = V * (jv(-n, S) - r ** 2 * C_n(n, S))
    g3 = V * (jv(-n, S) * (1 - 0.5 * W_n(n, S) * r ** 2) - r ** 2 * C_n(n, S))
    dstar = -np.pi * alpha * (alpha / aac) * B_n(n, S)
    return g1, g3C, g3, dstar


def exact(alpha, aac, S, n):
    """Minimum splitting / 2 of the two branches nearest 0, and its eta offset from n q."""
    q = 2 * np.pi * aac
    K = auto_K(S) + abs(n) + 6
    fl = FloquetLaue(alpha, aac, S, K)

    def split(d):
        w = np.linalg.eigvalsh(fl.H(np.array([n * q + d]))[0])
        two = np.sort(w[np.argsort(np.abs(w))[:2]])
        return two[1] - two[0]

    span = 3 * np.pi * alpha * (alpha / aac) * max(abs(B_n(n, S)), 0.5) + 0.5
    res = minimize_scalar(split, bounds=(-span, span), method="bounded", options={"xatol": 1e-10})
    return res.fun / 2, res.x


if __name__ == "__main__":
    print("== generic S: |g| exact vs 1st order, 3rd order (C only), 3rd order full; resonance shift ==")
    print(" n  S    R     | |g| exact    err g1      err g3(C)    err g3(full) | delta* exact   predicted")
    rows = []
    for n, S in ((0, 1.0), (1, 1.5), (1, 3.0), (2, 2.5), (2, 4.5)):
        for R in (0.01, 0.02, 0.05, 0.1, 0.2):
            aac = 10.0
            alpha = R * aac
            ge, de = exact(alpha, aac, S, n)
            g1, g3C, g3, ds = predictions(alpha, aac, S, n)
            rows.append((n, S, R, ge, abs(ge - abs(g1)), abs(ge - abs(g3C)), abs(ge - abs(g3)), de, ds))
            print(f" {n:+d} {S:3.1f} {R:5.2f} | {ge:.6e} {rows[-1][4]:.3e} {rows[-1][5]:.3e} {rows[-1][6]:.3e} |"
                  f" {de:+.6f} {ds:+.6f}", flush=True)
    print("\n== error scaling exponents (log-log slope over R = 0.01-0.1), relative to |g|/V ==")
    for n, S in ((0, 1.0), (1, 1.5), (1, 3.0), (2, 2.5), (2, 4.5)):
        sel = [r for r in rows if r[0] == n and r[1] == S and r[2] <= 0.1]
        R = np.array([r[2] for r in sel])
        V = np.pi * R * 10.0
        for idx, lab in ((4, "g1"), (5, "g3(C)"), (6, "g3(full)")):
            e = np.array([r[idx] for r in sel]) / V
            slope = np.polyfit(np.log(R), np.log(np.maximum(e, 1e-300)), 1)[0]
            print(f"  n={n} S={S}: {lab:9s} rel. error slope {slope:5.2f}  (values {np.array2string(e, precision=2)})")
        de = np.array([r[7] for r in sel]); ds = np.array([r[8] for r in sel])
        if np.any(np.abs(ds) > 0):
            sl = np.polyfit(np.log(R), np.log(np.maximum(np.abs(de - ds), 1e-300)), 1)[0]
            print(f"  n={n} S={S}: delta* residual slope {sl:5.2f}  (|delta*| slope "
                  f"{np.polyfit(np.log(R), np.log(np.maximum(np.abs(de), 1e-300)), 1)[0]:.2f})")
    print("\n== at Bessel zeros: |g| exact vs pi alpha (R/2)^2 |C_n|; renormalised zero dS vs R ==")
    for n in (0, 1, 2):
        j = jn_zeros(n, 1)[0]
        Cn = C_n(n, j)
        dJ = jvp(-n, j)
        for R in (0.02, 0.05, 0.1, 0.2):
            aac = 10.0
            alpha = R * aac
            ge, _ = exact(alpha, aac, j, n)
            gp = np.pi * alpha * (R / 2) ** 2 * abs(Cn)
            # renormalised zero: S where the exact splitting is minimal
            f = lambda S: exact(alpha, aac, S, n)[0]  # noqa: E731
            pred = (R / 2) ** 2 * Cn / dJ
            sp = max(3 * abs(pred), 1e-4)
            Sst = minimize_scalar(f, bounds=(j - sp, j + sp), method="bounded", options={"xatol": 1e-9}).x
            print(f"  n={n} R={R:5.2f}: |g| exact {ge:.4e} pred {gp:.4e} (ratio {ge/gp:.4f}) | "
                  f"dS exact {Sst-j:+.6e} pred {pred:+.6e} (ratio {(Sst-j)/pred:.4f}) | min |g| at S* {f(Sst):.2e}",
                  flush=True)
