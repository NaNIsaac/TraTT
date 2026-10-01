"""Exact harmonic (Floquet) solution of Laue two-beam TT for a uniform travelling grating.

For phi(zeta) = S sin(q zeta - q t) over the whole slab (crystal fully filled), expand

    D0 = sum_k a_k(zeta) e^{+i k q zeta},      Dh = sum_m b_m(zeta) e^{-i m q zeta}.

With e^{-i phi} = sum_p J_p(S) e^{-i p (q zeta - q t)} the TT equations become a
*constant-coefficient* system, x' = i H x, x = (a_k, b_m):

    H[a_k, a_k] = -k q                 H[b_m, b_m] = eta + m q
    H[a_k, b_m] = pi alpha J_{k+m}(S) e^{-i (k+m) q t}
    H[b_m, a_k] = conj(H[a_k, b_m])    (H is Hermitian)

The expansion is redundant at zeta = 0; the choice a_k(0) = delta_k0, b_m(0) = 0 satisfies
D0(0) = 1, Dh(0) = 0, so by ODE uniqueness the series solves the original problem exactly.
Truncation to |k|, |m| <= K is the only approximation (convergence is tested).

Channel eta = n q is the resonance of b_{-n} with a_0 (diagonals 0 and eta - n q), with
direct coupling J_{-n}; every other a_k, b_m is off resonance by a multiple of q.
The eigenvalues of H are the Floquet propagation constants (dispersion-surface branches).

`mask(k, m) -> bool` removes individual couplings a_k <-> b_m (path-removal diagnostics).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.special import jv

__all__ = ["FloquetLaue", "auto_K"]


def auto_K(S: float, eta_max_orders: float = 0.0, tol: float = 1e-10, extra: int = 3) -> int:
    """Smallest K with sum_{|p| > K} J_p(S)^2 < tol, widened to cover the eta window and a margin."""
    K = 0
    while 1.0 - (jv(0, S) ** 2 + 2 * sum(jv(p, S) ** 2 for p in range(1, K + 1))) > tol:
        K += 1
    return int(max(K, int(np.ceil(eta_max_orders))) + extra)


@dataclass
class FloquetLaue:
    alpha: float
    alpha_ac: float
    S: float
    K: int
    t: float = 0.0
    mask: Callable | None = None

    def __post_init__(self):
        self.q = 2.0 * np.pi * self.alpha_ac
        self.ks = np.arange(-self.K, self.K + 1)
        self.ms = np.arange(-self.K, self.K + 1)
        n = self.ks.size
        H0 = np.zeros((2 * n, 2 * n), dtype=np.complex128)
        H0[np.arange(n), np.arange(n)] = -self.ks * self.q
        for i, k in enumerate(self.ks):
            for j, m in enumerate(self.ms):
                if self.mask is not None and not self.mask(int(k), int(m)):
                    continue
                v = np.pi * self.alpha * jv(k + m, self.S) * np.exp(-1j * (k + m) * self.q * self.t)
                H0[i, n + j] = v
                H0[n + j, i] = np.conj(v)
        self.H0 = H0
        self.n = n

    def H(self, eta):
        """Stack of H(eta), shape (n_eta, 2n, 2n)."""
        eta = np.atleast_1d(np.asarray(eta, dtype=np.float64))
        H = np.broadcast_to(self.H0, (eta.size,) + self.H0.shape).copy()
        idx = self.n + np.arange(self.n)
        H[:, idx, idx] = eta[:, None] + self.ms[None, :] * self.q
        return H

    def propagate(self, eta, zeta=1.0, chunk=2048):
        """Harmonic amplitudes x(zeta) for every eta: shape (n_eta, 2n)."""
        eta = np.atleast_1d(np.asarray(eta, dtype=np.float64))
        out = np.empty((eta.size, 2 * self.n), dtype=np.complex128)
        x0 = np.zeros(2 * self.n, dtype=np.complex128)
        x0[self.K] = 1.0  # a_0
        for i in range(0, eta.size, chunk):
            w, U = np.linalg.eigh(self.H(eta[i:i + chunk]))
            c = np.conj(U[:, self.K, :])  # U^dagger x0
            out[i:i + chunk] = np.einsum("eij,ej->ei", U, np.exp(1j * w * zeta) * c)
        return out

    def fields(self, eta, zeta=1.0):
        """Physical exit amplitudes (D0, Dh) at depth zeta."""
        x = self.propagate(eta, zeta)
        a, b = x[:, : self.n], x[:, self.n:]
        D0 = a @ np.exp(1j * self.ks * self.q * zeta)
        Dh = b @ np.exp(-1j * self.ms * self.q * zeta)
        return D0, Dh

    def rocking(self, eta):
        return np.abs(self.fields(eta)[1]) ** 2

    def eigen(self, eta):
        """Floquet propagation constants and eigenvectors at each eta."""
        return np.linalg.eigh(self.H(eta))
