"""Entry 47b: map the dark loci of orders n = 0, 1, 2 in the (R, S) plane (alpha_ac = 2 and 5).

For each (alpha, S): phase-averaged (4 phases) window integrals p_n (|eta - n q| < q/2) and the
independent-channel values; ratio D_n = p_n / p_n^ind. Dark valleys D_n << 1 with |J_n| large are
multi-path dark states; valleys that start at Bessel zeros j_{n,k} as R -> 0 are renormalised zeros.
Third-order prediction of the renormalised zero: J_{-n}(S) = (R/2)^2 C_n(S).
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import sys  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.special import jv  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e43_phase import laue_line  # noqa: E402
from tratt.core import required_nz, solve_laue  # noqa: E402
from tratt.drive import Drive  # noqa: E402

NS = (0, 1, 2)


def point(args):
    alpha, aac, S = args
    q = 2 * np.pi * aac
    eta = np.arange(-2.5 * q, 2.5 * q, 0.15)
    nz = required_nz(alpha, aac, S=S, eta_max=2.5 * q)
    I = np.mean([np.abs(solve_laue(alpha, eta, phi=Drive(S=S, alpha_ac=aac).steady_phi(k / (4 * aac)),
                                   n_z=nz)[1]) ** 2 for k in range(4)], axis=0)
    out = []
    for n in NS:
        w = np.abs(eta - n * q) < q / 2
        ind = laue_line(eta - n * q, np.pi * alpha * abs(jv(n, S)))
        out.append(np.trapezoid(np.where(w, I, 0), eta) / max(np.trapezoid(np.where(w, ind, 0), eta), 1e-30))
    return out


if __name__ == "__main__":
    for aac in (2.0, 5.0):
        Rs = np.linspace(0.05, 1.5, 30)
        Ss = np.linspace(0.2, 6.0, 59)
        grid = [(R * aac, aac, S) for R in Rs for S in Ss]
        with Pool(11) as pool:
            D = np.array(pool.map(point, grid, chunksize=8)).reshape(len(Rs), len(Ss), len(NS))
        np.savez(HERE / f"e47b_darkmap_q{aac}.npz", Rs=Rs, Ss=Ss, D=D)
        print(f"aac={aac} done", flush=True)
        for j, n in enumerate(NS):
            print(f"  n={n}: darkest points with |J_n| > 0.3 (R, S, D):")
            cand = [(D[i, k, j], Rs[i], Ss[k]) for i in range(len(Rs)) for k in range(len(Ss)) if abs(jv(n, Ss[k])) > 0.3]
            for d, R, S in sorted(cand)[:6]:
                print(f"    R={R:.2f} S={S:.2f} D={d:.4f} |J_n|={abs(jv(n,S)):.3f}")
