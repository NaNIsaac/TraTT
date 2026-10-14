"""Entry 64: Sondhauss & Wark, Acta Cryst. A59, 7-13 (2003), eq. 29, against full TT and the Entry 59 form.

Single acoustic mode filling an absorbing symmetric-Bragg crystal (solve_bragg, alpha = 8, kappa = 6 as in
Entry 59), phi(zeta, t) = S cos(W (t - zeta)), quasi-static. Estimator: the Fourier component of R(t) at W,
    F(eta) = mean_t[(R - <R>) e^{iWt}] = conj(R0) dR_- + R0 conj(dR_+),
where dR_s is the e^{i s W t} part of the reflection amplitude at the surface, to first order in S.

Two-field form of the solve_bragg Riccati equation (a = D0, b = Dh, R = b/a, e = e^{i phi}):
    a' = c e b,   b' = -c conj(e) a - i eta_c b,   c = i pi alpha, eta_c = eta + i kappa.
Linearise about a perfect-crystal solution (a_p, b_p) with e ~ 1 + i phi, and take the mode ansatz
delta = (A_s, B_s) e^{i s W (t - zeta)} times a depth factor e^{lam zeta}; with L = lam - i s W:
    L A - c B           = i c (S/2) b_p
    c A + (L + i eta_c) B = i c (S/2) a_p
  SW   : S&W eq. 29 literally -- (a_p, b_p) = (1, R0) constant, lam = 0, satellite amplitude read as B
         (Dh = Dh^p + d_h, eq. 31). Their 2 Omega/omega terms are dropped (quasi-static, Entry 62).
  SWn  : same constant ansatz, but normalised to the incident wave, dR = B - R0 A.
  Bloch: (a_p, b_p) = (1, R0) e^{c R0 zeta} (the actual semi-infinite Bloch wave), lam = c R0,
         dR = B - R0 A. Must equal the Entry 59 Riccati form D_s = i c S (1 - R0^2) / (2 (Lambda - i s W)).
Usage: python e64_sw_compare.py
"""
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e59_transfer import root  # noqa: E402
from tratt.core import required_nz, solve_bragg  # noqa: E402

ALPHA, KAPPA, S = 8.0, 6.0, 0.02
C = 1j * np.pi * ALPHA
WD = 2 * np.pi * ALPHA  # Darwin half-width in eta


def modes(eta, W, kind):
    """dR_+, dR_- at the surface for one prediction kind ('SW', 'SWn', 'Bloch', 'Riccati')."""
    R0 = root(ALPHA, eta, KAPPA)
    eta_c = eta + 1j * KAPPA
    out = []
    for s in (+1, -1):
        if kind == "Riccati":
            lam_r = 1j * eta_c + 2 * C * R0
            out.append(1j * C * S * (1 - R0 ** 2) / (2 * (lam_r - 1j * s * W)))
            continue
        lam = C * R0 if kind == "Bloch" else 0.0
        L = lam - 1j * s * W
        rhs1, rhs2 = 1j * C * (S / 2) * R0, 1j * C * (S / 2)
        det = L * (L + 1j * eta_c) + C * C
        A = (rhs1 * (L + 1j * eta_c) + C * rhs2) / det
        B = (L * rhs2 - C * rhs1) / det
        out.append(B if kind == "SW" else B - R0 * A)
    return R0, out[0], out[1]


def predict(eta, W, kind):
    R0, dp, dm = modes(eta, W, kind)
    return np.conj(R0) * dm + R0 * np.conj(dp)


def tt(eta, W, nt=32, nz=None):
    t = 2 * np.pi * np.arange(nt) / (W * nt)  # one acoustic period, exact for the W component
    nz = nz or required_nz(ALPHA, W / (2 * np.pi), per_cycle=24, S=S, eta_max=np.abs(eta).max())

    def phi(zeta):
        return S * np.cos(W * (t[:, None] - zeta))
    R = np.abs(solve_bragg(ALPHA, eta, phi=phi, n_z=nz, kappa=KAPPA)) ** 2
    return np.mean((R - R.mean(axis=0)) * np.exp(1j * W * t[:, None]), axis=0), nz


def grid(W):
    g = [np.linspace(-2.5 * W, 2.5 * W, 201), np.linspace(-3 * WD, 3 * WD, 61)]
    g += [s * W + np.linspace(-3 * WD, 3 * WD, 61) for s in (+1, -1)]
    return np.unique(np.concatenate(g))


def stats(Ft, Fp):
    big = np.abs(Ft) > 0.02 * np.abs(Ft).max()
    r = np.abs(Fp[big]) / np.abs(Ft[big])
    cerr = np.abs(Fp[big] - Ft[big]) / np.abs(Ft[big])  # complex (phase-sensitive) error
    return np.median(r), r.min(), r.max(), np.median(cerr), cerr.max()


if __name__ == "__main__":
    # analytic identity first: Bloch-consistent two-field solution == Entry 59 Riccati form
    e = np.linspace(-300, 300, 41)
    for W in (WD / 0.05, WD / 1.0):
        ref = predict(e, W, "Riccati")  # F vanishes at eta = 0 by symmetry: compare to max |F|
        d = np.max(np.abs(predict(e, W, "Bloch") - ref)) / np.abs(ref).max()
        print(f"identity check W={W:.0f}: max |Bloch - Riccati| / max|F| = {d:.1e}")
    out = {}
    for RD in (0.05, 0.2, 0.5, 1.0, 2.0):
        W = WD / RD
        eta = grid(W)
        Ft, nz = tt(eta, W)
        print(f"\nR_D = {RD} (W = {W:.0f} = {1 / RD:.1f} Darwin half-widths), n_z = {nz}")
        if RD == 0.5:
            Ft2, _ = tt(eta, W, nt=64, nz=2 * nz)
            print(f"  convergence (nt x2, n_z x2): max |dF| / max|F| = "
                  f"{np.max(np.abs(Ft2 - Ft)) / np.abs(Ft).max():.1e}")
        rec = {"eta": eta.tolist(), "W": W, "TT": [Ft.real.tolist(), Ft.imag.tolist()]}
        for kind in ("SW", "SWn", "Riccati"):
            Fp = predict(eta, W, kind)
            med, lo, hi, cmed, cmax = stats(Ft, Fp)
            print(f"  {kind:8s} |F_pred|/|F_TT| median {med:.3f} range [{lo:.3f}, {hi:.3f}];"
                  f"  complex rel err median {cmed:.3f} max {cmax:.2f}")
            rec[kind] = [Fp.real.tolist(), Fp.imag.tolist()]
        # the point a measurer reads: the satellite, eta = +-W
        for x in (W, -W):  # F(eta = 0) = 0 by symmetry, so no ratio there
            i = np.argmin(np.abs(eta - x))
            print(f"   eta = {eta[i]:+8.1f}: |F_TT| {abs(Ft[i]):.3e}  SW/TT {abs(predict(eta[i:i+1], W, 'SW')[0]) / abs(Ft[i]):.3f}"
                  f"  SWn/TT {abs(predict(eta[i:i+1], W, 'SWn')[0]) / abs(Ft[i]):.3f}"
                  f"  Riccati/TT {abs(predict(eta[i:i+1], W, 'Riccati')[0]) / abs(Ft[i]):.3f}")
        out[f"RD{RD}"] = rec
    json.dump(out, open(HERE / "e64_sw.json", "w"))
