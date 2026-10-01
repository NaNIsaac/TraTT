"""Entry 45: Bessel-zero breakdown — can a nominally dark Floquet order come back?

Uniform travelling wave, filled crystal, t = 0. At S = j (a zero of J_n), the direct coupling
of channel eta = n q (a_0 <-> b_{-n}) vanishes. H couples only a's to b's (bipartite), so the
lowest indirect path is third order, a_0 -> b_m -> a_k -> b_{-n}. Lowdin / adiabatic
elimination of the off-resonant a_k (energy -k q) and b_m (energy (m + n) q) gives

    g_eff = pi alpha [ J_{-n}(S) - (pi alpha / q)^2 C_n(S) ],
    C_n(S) = sum_{k != 0, m != -n} J_{k-n}(S) J_{k+m}(S) J_m(S) / (k (m + n)),

so at the zero the channel population is sin^2|g_eff| ~ (pi alpha)^2 (alpha / 2 alpha_ac)^4 C_n^2.
Measured: max over eta near n q of the Floquet component |b_{-n}(1)|^2 (resonance shifted by
second-order level shifts), the physical window integral and window maximum, for integer and
half-integer alpha_ac; plus path removal of the dominant (m, k) terms.
Usage: python e45_zero.py [scan] [paths] [phys]
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.special import jn_zeros, jv  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402

ZEROS = {0: jn_zeros(0, 1)[0], 1: jn_zeros(1, 1)[0], 2: jn_zeros(2, 1)[0]}


def C_n(n, S, K=40):
    tot = 0.0
    for k in range(-K, K + 1):
        if k == 0:
            continue
        for m in range(-K, K + 1):
            if m == -n:
                continue
            tot += jv(k - n, S) * jv(k + m, S) * jv(m, S) / (k * (m + n))
    return tot


def terms(n, S, K=12):
    out = []
    for k in range(-K, K + 1):
        for m in range(-K, K + 1):
            if k == 0 or m == -n:
                continue
            out.append((jv(k - n, S) * jv(k + m, S) * jv(m, S) / (k * (m + n)), m, k))
    return sorted(out, key=lambda x: -abs(x[0]))


def channel_population(alpha, aac, S, n, mask=None, width=None):
    q = 2 * np.pi * aac
    K = auto_K(S) + abs(n) + 3 + int(np.ceil(4 * np.pi * alpha / q))
    fl = FloquetLaue(alpha, aac, S, K, mask=mask)
    w = width or min(0.45 * q, 12 + 3 * (np.pi * alpha) ** 2 / q)
    eta = n * q + np.linspace(-w, w, 1601)
    x = fl.propagate(eta)
    b = x[:, fl.n + (fl.K - n)]  # component b_{-n}
    P = np.abs(b) ** 2
    k = int(np.argmax(P))
    return P[k], eta[k] - n * q, fl, eta


def section_scan():
    print("== channel population at the Bessel zero vs third-order prediction ==")
    print(" n  S      alpha  aac   | max|b_-n|^2   pred sin^2|g_eff|   ratio | resonance shift")
    for n, S in ZEROS.items():
        Cn = C_n(n, S)
        print(f"  C_{n}(S={S:.4f}) = {Cn:+.5f}")
        for alpha in (0.5, 1.0, 2.0):
            for aac in (40.0, 20.0, 10.0, 5.0, 3.0, 2.0, 1.5, 1.0):
                P, sh, _, _ = channel_population(alpha, aac, S, n)
                g = np.pi * alpha * (jv(-n, S) - (alpha / (2 * aac)) ** 2 * Cn)
                pred = np.sin(abs(g)) ** 2
                print(f" {n:+d} {S:.3f} {alpha:5.2f} {aac:5.1f} | {P:12.3e} {pred:12.3e} {P / pred:9.3f} | {sh:+.3f}")


def section_paths():
    print("== which (m, k) paths carry the resurrection: largest terms of C_n and removal tests ==")
    for n, S in ZEROS.items():
        T = terms(n, S)
        print(f" n={n} S={S:.4f}: top terms (value, m, k): " +
              ", ".join(f"({v:+.4f}, m={m}, k={k})" for v, m, k in T[:6]))
        alpha, aac = 1.0, 3.0
        P0, _, _, _ = channel_population(alpha, aac, S, n)
        # remove the dominant intermediate b_m entirely (all its couplings) and the dominant a_k
        _, m1, k1 = T[0]
        Pm, _, _, _ = channel_population(alpha, aac, S, n, mask=lambda k, m: m != m1)
        Pk, _, _, _ = channel_population(alpha, aac, S, n, mask=lambda k, m: k != k1)
        # keep only the dominant path's three couplings + the resonant pair
        keep = {(0, m1), (k1, m1), (k1, -n), (0, -n)}
        Pkeep, _, _, _ = channel_population(alpha, aac, S, n, mask=lambda k, m: (k, m) in keep)
        g1 = np.pi * alpha * (alpha / (2 * aac)) ** 2 * T[0][0]
        print(f"   alpha={alpha} aac={aac}: full {P0:.3e} | without b_{m1}: {Pm:.3e} | without a_{k1}: {Pk:.3e}"
              f" | only path (m={m1},k={k1}): {Pkeep:.3e} (pred sin^2 of that term {np.sin(abs(g1))**2:.3e})")


def section_phys():
    print("== physical intensity near eta = n q at the zero: window max, window integral, Floquet component ==")
    for n, S in ((1, ZEROS[1]),):
        for alpha in (1.0, 2.0):
            for aac in (20.0, 20.5, 5.0, 5.5, 2.0, 2.5, 1.0):
                q = 2 * np.pi * aac
                P, sh, fl, eta = channel_population(alpha, aac, S, n, width=0.45 * q)
                D0, Dh = fl.fields(eta)
                I = np.abs(Dh) ** 2
                ref = FloquetLaue(alpha, aac, 1.0, fl.K)  # same alpha, a bright S for scale
                Ir = np.abs(ref.fields(eta)[1]) ** 2
                print(f"  n={n} alpha={alpha} aac={aac:5.1f}: |b_-n|^2 max {P:.2e} | phys window max {I.max():.2e}"
                      f" int {np.trapezoid(I, eta):.3e} | same window at S=1: max {Ir.max():.2f} int {np.trapezoid(Ir, eta):.2f}")


if __name__ == "__main__":
    for s in (sys.argv[1:] or ["scan", "paths", "phys"]):
        globals()[f"section_{s}"]()
