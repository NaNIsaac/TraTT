"""Entry 46: exact two-band Floquet structure ("X-ray band structure of the acoustic superlattice").

The TT system is a 2x2 ODE with coefficients periodic in zeta with period Lam = 2 pi / q.
Its monodromy M(eta) over one period has eigenvalues e^{i lambda_{1,2} Lam}; lambda_1 + lambda_2 = eta
(mod q) and the pair is unitary. Reduced quasi-propagation constants lambda in (-q/2, q/2].
  gap(eta)   = min distance between the two branches (mod q)
  Pendelloesung over the slab is set by Delta lambda = lambda_1 - lambda_2.
Independent-channel picture: gaps of width 2 pi alpha |J_n| at eta = n q, isolated.
Also: eigenvector mixing (how much each Bloch mode is 'forward' vs 'diffracted').
Outputs: e46_bands.npz for the figure; printed gap table.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.special import jv  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.core import solve_laue  # noqa: E402
from tratt.drive import Drive  # noqa: E402


def monodromy(alpha, aac, S, eta, nz=4000):
    """M(eta) over one acoustic period (in units of the slab thickness, Lam = 1 / aac)."""
    lam = 1.0 / aac
    d = Drive(S=S, alpha_ac=aac)
    phi = d.steady_phi(0.0)
    # integrate the slab problem on [0, lam] by rescaling: zeta = lam * u, u in [0, 1]
    # equivalent coupling alpha * lam and detuning eta * lam with phase phi(lam u)
    from tratt import core
    cols = []
    for D0i, Dhi in ((1.0, 0.0), (0.0, 1.0)):
        D0 = np.full(eta.shape, D0i, dtype=complex)
        Dh = np.full(eta.shape, Dhi, dtype=complex)
        c = 1j * np.pi * alpha * lam
        h = 1.0 / nz

        def f(u, d0, dh):
            e = np.exp(1j * phi(np.float64(lam * u)))
            return c * e * dh, c * np.conj(e) * d0 + 1j * eta * lam * dh

        for i in range(nz):
            u = i * h
            k1 = f(u, D0, Dh)
            k2 = f(u + h / 2, D0 + h / 2 * k1[0], Dh + h / 2 * k1[1])
            k3 = f(u + h / 2, D0 + h / 2 * k2[0], Dh + h / 2 * k2[1])
            k4 = f(u + h, D0 + h * k3[0], Dh + h * k3[1])
            D0 = D0 + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            Dh = Dh + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        cols.append(np.stack([D0, Dh], -1))
    return np.stack(cols, -1)  # (n_eta, 2, 2)


def bands(alpha, aac, S, eta, nz=None):
    q = 2 * np.pi * aac
    nz = nz or int(60 + 40 * (S + np.abs(eta).max() / q + alpha / aac))
    M = monodromy(alpha, aac, S, eta, nz)
    mu, V = np.linalg.eig(M)
    lam = np.angle(mu) * aac  # quasi-propagation constant, (-q/2, q/2]
    order = np.argsort(lam, axis=1)
    lam = np.take_along_axis(lam, order, 1)
    V = np.take_along_axis(V, order[:, None, :], 2)
    dl = lam[:, 1] - lam[:, 0]
    gap = np.minimum(dl, q - dl)
    mix = np.abs(V[:, 1, :]) ** 2 / np.sum(np.abs(V) ** 2, axis=1)  # diffracted weight of each mode
    return lam, gap, mix, np.abs(np.linalg.det(M))


if __name__ == "__main__":
    out = {}
    print(" S    aac  alpha   R   | gaps at eta = n q (n = 0..3): measured / 2 pi alpha |J_n| | max gap / q | |det M|-1")
    for S in (1.5, 3.8317):
        for aac, alpha in ((5.0, 0.25), (5.0, 1.0), (5.0, 2.5), (5.0, 5.0), (2.0, 4.0)):
            q = 2 * np.pi * aac
            eta = np.linspace(-0.5 * q, 3.5 * q, 2001)
            lam, gap, mix, det = bands(alpha, aac, S, eta)
            gs = [gap[np.argmin(np.abs(eta - n * q))] for n in range(4)]
            pred = [2 * np.pi * alpha * abs(jv(n, S)) for n in range(4)]
            print(f" {S:4.2f} {aac:4.1f} {alpha:5.2f} {alpha/aac:4.2f} | " +
                  " ".join(f"{g:6.2f}/{p:6.2f}" for g, p in zip(gs, pred)) +
                  f" | {gap.max()/q:.3f} | {np.abs(det-1).max():.1e}", flush=True)
            out[f"S{S}_q{aac}_a{alpha}"] = np.array([eta, lam[:, 0], lam[:, 1], gap, mix[:, 0]])
    np.savez(HERE / "e46_bands.npz", **out)
