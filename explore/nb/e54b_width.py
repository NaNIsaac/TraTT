"""Entry 54b: target contribution vs integration half-width w and period mismatch eps.

alpha = 1, S = j_11, channel n = 1 (target b_-1), alpha_ac = 20 + eps. Exact phase-averaged intensity
<I(eta)> = sum_m |b_m(eta)|^2 (Entry 53). For a centre c (nominal eta = q, or the numerical resonance
q + delta*, delta* = location of the minimum two-branch gap, as in e52_heff.exact), integrate the
target |b_-1|^2 and the total over [c - w, c + w] on a grid that is dense near c, for w from 1e-8 q
to 0.45 q, and report the fraction F(w) = int target / int total and the w at which F = 0.5, 0.9.
Grid convergence: the same with twice the points.
Usage: python e54b_width.py [base] [eps] [map]
"""
import os
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1" if "map" in sys.argv else "8")
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.integrate import cumulative_trapezoid  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402
from scipy.special import jn_zeros  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402

J11 = jn_zeros(1, 1)[0]
ALPHA, N = 1.0, 1


def resonance(fl, q):
    def split(d):
        w = np.linalg.eigvalsh(fl.H(np.array([N * q + d]))[0])
        two = np.sort(w[np.argsort(np.abs(w))[:2]])
        return (two[1] - two[0]) / 2
    return minimize_scalar(split, bounds=(-0.05, 0.05), method="bounded", options={"xatol": 1e-12}).x


def profile(aac, center, npts=1500, K=None):
    q = 2 * np.pi * aac
    K = K or auto_K(J11) + 8
    fl = FloquetLaue(ALPHA, aac, J11, K)
    u = np.concatenate([[0.0], np.geomspace(1e-8 * q, 0.45 * q, npts)])
    eta = np.concatenate([center - u[::-1][:-1], center + u])
    x = fl.propagate(eta)
    b = x[:, fl.n:]
    tgt = np.abs(b[:, list(fl.ms).index(-N)]) ** 2
    tot = np.sum(np.abs(b) ** 2, axis=1)
    i0 = u.size - 1
    # integrals over [c - w, c + w] for each w = u
    def sym(y):
        right = cumulative_trapezoid(y[i0:], eta[i0:], initial=0.0)
        left = cumulative_trapezoid(y[i0::-1][::1], -eta[i0::-1], initial=0.0)
        return right + left
    return u, sym(tgt), sym(tot), fl, q


def crossings(u, F, q):
    out = {}
    for thr in (0.5, 0.9):
        bad = np.flatnonzero(~(F[1:] >= thr)) + 1       # F >= thr contiguously from the centre up to w
        k = bad[0] - 1 if bad.size else u.size - 1
        out[thr] = float(u[k] / q) if k >= 1 else None
    return out


def section_base():
    aac = 20.0
    q = 2 * np.pi * aac
    fl = FloquetLaue(ALPHA, aac, J11, auto_K(J11) + 8)
    dstar = resonance(fl, q)
    print(f"alpha_ac = 20: numerical resonance delta* = {dstar:+.4e}")
    res = {}
    for lab, c in (("nominal", N * q), ("resonance", N * q + dstar)):
        for npts, KK in ((1500, None), (3000, None), (1500, auto_K(J11) + 12)):
            u, It, Ia, _, _ = profile(aac, c, npts, KK)
            F = np.where(Ia > 0, It / Ia, np.nan)
            cr = crossings(u, F, q)
            res[(lab, npts)] = (u, It, Ia, F)
            sel = [np.argmin(abs(u - w * q)) for w in (1e-6, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.45)]
            print(f" {lab:9s} npts {npts} K {KK}: F at w/q = " + " ".join(f"{u[i]/q:.0e}:{F[i]:.4f}" for i in sel) +
                  f" | w/q at F=0.5: {cr[0.5]}, F=0.9: {cr[0.9]}")
            for thr in (0.5, 0.9):
                if cr[thr]:
                    i = np.argmin(abs(u - cr[thr] * q))
                    print(f"     F >= {thr}: w/q <= {cr[thr]:.4e}; there int target {It[i]:.3e}, int total {Ia[i]:.3e}")
    json.dump({f"{k[0]}_{k[1]}": dict(wq=(v[0] / q).tolist(), It=v[1].tolist(), Ia=v[2].tolist())
               for k, v in res.items() if k[1] == 3000}, open(HERE / "e54b_base.json", "w"))


def section_eps():
    eps = np.concatenate([-np.geomspace(0.05, 1e-4, 13), [0.0], np.geomspace(1e-4, 0.05, 13)])
    rows = []
    for e in eps:
        aac = 20.0 + e
        q = 2 * np.pi * aac
        fl = FloquetLaue(ALPHA, aac, J11, auto_K(J11) + 8)
        dstar = resonance(fl, q)
        u, It, Ia, _, _ = profile(aac, N * q + dstar, 1500)
        F = np.where(Ia > 0, It / Ia, np.nan)
        cr = crossings(u, F, q)
        rows.append(dict(eps=float(e), dstar=float(dstar), wq=(u / q).tolist(), F=np.nan_to_num(F).tolist(),
                         It=It.tolist(), Ia=Ia.tolist(), w50=cr[0.5], w90=cr[0.9]))
        print(f" eps={e:+.4e}: delta*={dstar:+.3e} | w/q (F>=0.5) {cr[0.5]} | (F>=0.9) {cr[0.9]} | "
              f"F at w=1e-3 q {F[np.argmin(abs(u - 1e-3 * q))]:.3f}, at 0.45 q {F[-1]:.2e}", flush=True)
    json.dump(rows, open(HERE / "e54b_eps.json", "w"))


def section_res():
    """Minimum-gap position vs eta = nq: truncation (K) and search tolerance (xatol) checks."""
    for aac in (20.0, 20.25, 10.0):
        q = 2 * np.pi * aac
        for K in (auto_K(J11) + 4, auto_K(J11) + 8, auto_K(J11) + 12):
            fl = FloquetLaue(ALPHA, aac, J11, K)

            def split(d):
                w = np.linalg.eigvalsh(fl.H(np.array([N * q + d]))[0])
                two = np.sort(w[np.argsort(np.abs(w))[:2]])
                return (two[1] - two[0]) / 2
            for tol in (1e-9, 1e-12):
                d = minimize_scalar(split, bounds=(-0.05, 0.05), method="bounded", options={"xatol": tol}).x
                print(f" alpha_ac {aac:6.2f} K {K} xatol {tol:.0e}: delta* = {d:+.4e}, g(delta*) = {split(d):.10e}, "
                      f"g(0) - g(delta*) = {split(0.0) - split(d):.2e}, delta*/q = {d / q:+.1e}")


def map_point(e):
    aac = 20.0 + e
    q = 2 * np.pi * aac
    dstar = resonance(FloquetLaue(ALPHA, aac, J11, auto_K(J11) + 8), q)
    u, It, Ia, _, _ = profile(aac, N * q + dstar, 1500)
    with np.errstate(invalid="ignore", divide="ignore"):
        F = It / Ia
    cr = crossings(u, F, q)
    return dict(eps=float(e), dstar=float(dstar), F=np.nan_to_num(F).tolist(), Ia=Ia.tolist(),
                w50=cr[0.5], w90=cr[0.9])


def section_map():
    """Linear eps grid for the summary figure; u/q is the same grid for every alpha_ac."""
    from multiprocessing import Pool
    eps = np.round(np.arange(-0.008, 0.00801, 2.5e-4), 6)
    with Pool(11) as pool:
        rows = pool.map(map_point, eps)
    q = 2 * np.pi * 20.0
    wq = (np.concatenate([[0.0], np.geomspace(1e-8 * q, 0.45 * q, 1500)]) / q).tolist()
    for r in rows:
        print(f" eps={r['eps']:+.5f}: w/q(F>=0.5) {r['w50']}, w/q(F>=0.9) {r['w90']}", flush=True)
    json.dump(dict(wq=wq, rows=rows), open(HERE / "e54b_map.json", "w"))


if __name__ == "__main__":
    for s in (sys.argv[1:] or ["base", "eps"]):
        globals()[f"section_{s}"]()
