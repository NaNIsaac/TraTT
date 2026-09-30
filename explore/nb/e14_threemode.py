"""Entry 14: a closed three-mode, two-layer model for I(eta, t) near satellite n.

Laue, step onset, t <= 1. e^{-i phi} = sum_m J_m(S) e^{-i m (q zeta - q t)}; channel
eta = n q is driven by harmonic m = -n. Keep only the harmonics that matter near
eta = n q: m = 0 (the main-line tail) and m = -n (the satellite).

Filled layer [0, t], with Dh = H0 + Hn e^{i n q zeta}, rotating-wave:
    D0' = i pi alpha (J_0 H0 + J_{-n} e^{+i n q t} Hn)
    H0' = i pi alpha J_0 D0            + i eta H0
    Hn' = i pi alpha J_{-n} e^{-i n q t} D0 + i (eta - n q) Hn
constant coefficients -> one matrix exponential.
Unfilled layer [t, 1]: unmodulated Laue, 2x2, exact.
Output Dh(1) = [layer 2] (D0(t), H0(t) + Hn(t) e^{i n q t}).

Compared with full TT on the window |eta - n q| < q/2.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402


def model(alpha, S, q, n, eta, t, keep_tail=True):
    c = 1j * np.pi * alpha
    J0, Jn = jv(0, S), jv(-n, S)
    ph = np.exp(1j * n * q * t)
    out = np.empty(eta.size, dtype=complex)
    for k, e in enumerate(eta):
        j0 = J0 if keep_tail else 0.0
        M1 = np.array([[0, c * j0, c * Jn * ph],
                       [c * j0, 1j * e, 0],
                       [c * Jn / ph, 0, 1j * (e - n * q)]])
        D0, H0, Hn = expm(M1 * t) @ np.array([1, 0, 0], dtype=complex)
        M2 = np.array([[0, c], [c, 1j * e]])
        out[k] = (expm(M2 * (1 - t)) @ np.array([D0, H0 + Hn * ph]))[1]
    return np.abs(out) ** 2


if __name__ == "__main__":
    for S, alpha, n in [(1.5, 2.0, 1), (0.5, 2.0, 1), (1.5, 1.0, -1), (3.0, 2.0, 2)]:
        d = Drive(S=S, alpha_ac=20.0)
        q = d.q
        eta = np.linspace((n - 0.5) * q, (n + 0.5) * q, 801)
        t = np.array([0.15, 0.3137, 0.5, 0.7037, 0.9])
        I = rocking_cube(alpha, d, eta, t, n_z=1500)
        row = []
        for k, tt in enumerate(t):
            Im = model(alpha, S, q, n, eta, tt)
            In = model(alpha, S, q, n, eta, tt, keep_tail=False)
            nrm = I[k].max()
            row.append(f"t={tt:.2f}: 3-mode {np.max(np.abs(Im - I[k]))/nrm:.3f} "
                       f"no-tail {np.max(np.abs(In - I[k]))/nrm:.3f}")
        print(f"S={S} alpha={alpha} n={n}  max|model-TT|/peak\n   " + "\n   ".join(row))
