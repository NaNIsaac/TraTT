"""Entry 49: transient formation (travelling step front) in five regimes of the phase diagram.

alpha_ac = 5. Window integrals p_n(t) for n = 0, 1, 2 (|eta - nq| < q/2), averaged over each acoustic
period (8 samples per period) to remove micromotion; compared with the independent-channel
formation law p_n^ind(t) = window integral of the Laue line of thickness t, coupling pi alpha |J_n|.
Fourier analysis of p_n(t) for t < 1: dominant frequencies (cycles per transit) vs the independent
expectation ~ alpha |J_n|.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.special import jn_zeros, jv  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.core import required_nz  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

AAC = 5.0
CASES = {
    "1 resolved (R 0.1)": (0.5, 1.5),
    "2 onset (R 0.5)": (2.5, 1.5),
    "3 strong overlap (R 1)": (5.0, 1.5),
    "4 J1 = 0, resurrected n=1 (R 1)": (5.0, float(jn_zeros(1, 1)[0])),
    "5 dark main line (R 1.05, S 0.8)": (5.25, 0.8),
}


def laue_thick(d, a, t):
    s = np.sqrt(a * a + d * d / 4)
    return a * a * np.sin(s * t) ** 2 / np.where(s > 0, s * s, 1)


def peaks_of(f, P, k=5, floor=0.05):
    i = [j for j in range(1, P.size - 1) if P[j] > P[j - 1] and P[j] >= P[j + 1] and P[j] > floor * P.max()]
    return [float(f[j]) for j in sorted(i, key=lambda j: -P[j])[:k]]


if __name__ == "__main__":
    # NOTE: a first version averaged over each acoustic period (8 samples); the period is 0.2 t_tr,
    # which coarse-grained time to 0.2 and quantised every frequency. Kept out; see notebook.
    q = 2 * np.pi * AAC
    out = {}
    for name, (alpha, S) in CASES.items():
        eta = np.arange(-2.5 * q, 2.5 * q, 0.12)
        t = np.arange(0, 1.0, 1 / (16 * AAC)) + 0.0037
        nz = required_nz(alpha, AAC, S=S, eta_max=2.5 * q)
        I = rocking_cube(alpha, Drive(S=S, alpha_ac=AAC), eta, t, n_z=nz)
        res = []
        for n in (0, 1, 2):
            w = np.abs(eta - n * q) < q / 2
            p = np.trapezoid(np.where(w, I, 0), eta, axis=1)
            c = I[:, np.argmin(np.abs(eta - n * q))]
            a = np.pi * alpha * abs(jv(n, S))
            pind = np.array([np.trapezoid(np.where(w, laue_thick(eta - n * q, a, x), 0), eta) for x in t])
            f = np.fft.rfftfreq(t.size, t[1] - t[0])
            Pp = np.abs(np.fft.rfft((p - p.mean()) * np.hanning(t.size)))
            Pc = np.abs(np.fft.rfft((c - c.mean()) * np.hanning(t.size)))
            res.append((n, p, pind, c))
            err = np.trapezoid(np.abs(p - pind), t) / max(np.trapezoid(pind, t), 1e-30)
            print(f"{name:34s} n={n}: alpha|J_n|={alpha*abs(jv(n,S)):.2f} | dev from independent formation {err:.3f} | "
                  f"integral spectrum peaks {np.round(peaks_of(f, Pp), 1).tolist()} | centre spectrum peaks "
                  f"{np.round(peaks_of(f, Pc), 1).tolist()}", flush=True)
        out[name] = (t, eta, res)
    np.save(HERE / "e49_time.npy", np.array(out, dtype=object), allow_pickle=True)
