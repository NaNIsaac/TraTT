"""Entry 44: what fails first when Floquet channels overlap? Reduced coupled-mode models vs full TT.

Models (Floquet basis a_k: forward-beam harmonics, b_m: diffracted harmonics):
  IND   independent channels: each b_{-n} with a_0 alone (two-level), intensities added
  STAR  a_0 + all b_m          (one shared forward beam, no forward sidebands)
  A1    a_{0,+-1} + all b_m
  A2    a_{0,+-1,+-2} + all b_m
  FULL  |k| <= K converged (= full TT to 1e-5)
Metrics against FULL: E_pop (order-window populations), E_spec (whole line), total ratio.
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
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402
from e43_phase import laue_line  # noqa: E402


def windows(eta, I, q, ns):
    return np.array([np.trapezoid(np.where(np.abs(eta - n * q) < q / 2, I, 0), eta) for n in ns])


def run(alpha, aac, S, deta=0.25):
    q = 2 * np.pi * aac
    N = int(np.ceil(S)) + 4
    half = N * q + 4 * np.pi * alpha + 20
    eta = np.arange(-half, half, deta)
    K = max(auto_K(S), N) + 3 + int(np.ceil(4 * np.pi * alpha / q))
    ns = np.arange(-N, N + 1)
    full = FloquetLaue(alpha, aac, S, K).rocking(eta)
    out = {"FULL": full}
    out["IND"] = sum(laue_line(eta - n * q, np.pi * alpha * abs(jv(n, S))) for n in ns)
    for name, kmax in (("STAR", 0), ("A1", 1), ("A2", 2)):
        out[name] = FloquetLaue(alpha, aac, S, K, mask=lambda k, m, kmax=kmax: abs(k) <= kmax).rocking(eta)
    pf = windows(eta, full, q, ns)
    res = {}
    for name, I in out.items():
        p = windows(eta, I, q, ns)
        res[name] = (np.abs(p - pf).sum() / pf.sum(), np.trapezoid(np.abs(I - full), eta) / np.trapezoid(full, eta),
                     np.trapezoid(I, eta) / np.trapezoid(full, eta))
    return res


if __name__ == "__main__":
    print(" aac  alpha   R  S  | model: E_pop / E_spec / total ratio (vs FULL)")
    for S in (1.5, 4.0):
        for aac in (5.0, 2.0):
            for R in (0.05, 0.2, 0.5, 1.0, 2.0):
                alpha = R * aac
                r = run(alpha, aac, S)
                print(f" {aac:4.1f} {alpha:5.2f} {R:4.2f} {S:3.1f} | " +
                      " | ".join(f"{k}: {v[0]:.3f}/{v[1]:.3f}/{v[2]:.2f}" for k, v in r.items() if k != "FULL"),
                      flush=True)
