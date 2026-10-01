"""Entry 48: spectral-peak trajectories vs R, and the exact Bloch-mode mixing points.

alpha_ac = 2, S in {0.5, 1.0, 2.3}, R = alpha / alpha_ac from 0.05 to 2. For each R:
  - phase-averaged (8 phases) full-TT spectrum over |eta| < 3.5 q; maxima > 5 % of max, in units of q
  - two-band Floquet structure (monodromy): eta where the diffracted weight of a Bloch mode
    crosses 1/2 (= centre of an avoided crossing of the dressed branches)
Distinguishes "dark because the coupling cancels" from "dark because the line moved".
Saves e48_traj_S{S}.npz.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e46_bands import bands  # noqa: E402
from tratt.core import required_nz, solve_laue  # noqa: E402
from tratt.drive import Drive  # noqa: E402

AAC = 2.0


def one(args):
    R, S = args
    alpha = R * AAC
    q = 2 * np.pi * AAC
    eta = np.arange(-3.5 * q, 3.5 * q, 0.08)
    nz = required_nz(alpha, AAC, S=S, eta_max=3.5 * q)
    I = np.mean([np.abs(solve_laue(alpha, eta, phi=Drive(S=S, alpha_ac=AAC).steady_phi(k / (8 * AAC)),
                                   n_z=nz)[1]) ** 2 for k in range(8)], axis=0)
    pk = [i for i in range(1, I.size - 1) if I[i] > I[i - 1] and I[i] >= I[i + 1] and I[i] > 0.05 * I.max()]
    e = np.linspace(0.0, 0.5 * q, 201)
    _, _, mix, _ = bands(alpha, AAC, S, e)
    m = mix[:, 0]
    cross = [float(e[i] / q) for i in range(1, e.size) if (m[i] - 0.5) * (m[i - 1] - 0.5) < 0]
    return R, [float(eta[i] / q) for i in pk], [float(I[i]) for i in pk], cross, I[::5]


if __name__ == "__main__":
    for S in (0.5, 1.0, 2.3):
        Rs = np.linspace(0.05, 2.0, 40)
        with Pool(11) as pool:
            res = pool.map(one, [(R, S) for R in Rs])
        np.savez(HERE / f"e48_traj_S{S}.npz", Rs=Rs, peaks=np.array([r[1] for r in res], dtype=object),
                 heights=np.array([r[2] for r in res], dtype=object), cross=np.array([r[3] for r in res], dtype=object),
                 spec=np.array([r[4] for r in res]), allow_pickle=True)
        print(f"S={S}:")
        for R, pk, h, cr, _ in res[::4]:
            inner = sorted(p for p in pk if 0 <= p <= 1.6)
            print(f"  R={R:.2f}: peaks in [0,1.6] q: {np.round(inner, 3).tolist()} | Bloch mixing points in [0, q/2]: "
                  f"{np.round(cr, 3).tolist()}", flush=True)
