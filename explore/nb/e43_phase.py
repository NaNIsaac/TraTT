"""Entry 43: phase diagram of the resolved-Floquet-channel picture (uniform travelling wave,
crystal fully filled, acoustic phase t = 0).

Independent-channel prediction: every order n is a separate Laue slab of thickness 1 and
coupling a_n = pi alpha |J_n(S)| centred at eta = n q, the orders adding *in intensity*:
    I_ind(eta) = sum_n a_n^2 sin^2(sqrt(a_n^2 + d^2/4)) / (a_n^2 + d^2/4),  d = eta - n q.
Metrics per point (cached in cache_e43/):
  E_spec = int |I - I_ind| / int I_ind              (whole-spectrum error)
  E_pop  = sum_n |p_n - p_n^ind| / sum_n p_n^ind     (order-population error, windows |eta-nq|<q/2)
  E_ctr  = max_n |I(nq) - sin^2(a_n)| over orders with |J_n| > 0.02
  PR, PR_ind   participation ratio over order windows; tot, tot_ind  total diffracted intensity
  npk    spectral maxima above 1 % of the maximum
Solver: full TT (RK4) or the exact Floquet expansion, whichever is cheaper (both validated
to ~1e-5, tests/test_floquet.py).
Usage: python e43_phase.py            (full scan)   |  python e43_phase.py point a aac S
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json  # noqa: E402
import sys  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.special import jv  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.core import required_nz, solve_laue  # noqa: E402
from tratt.drive import Drive  # noqa: E402
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402

CACHE = HERE / "cache_e43"
ALPHAS = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
AACS = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 40.0]
SS = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0]


def laue_line(d, a):
    s2 = a * a + d * d / 4
    return np.where(s2 > 0, a * a * np.sin(np.sqrt(s2)) ** 2 / np.where(s2 > 0, s2, 1), 0.0)


def spectrum(alpha, aac, S, deta=0.3, t=0.0, solver="auto", Kextra=0):
    q = 2 * np.pi * aac
    N = int(np.ceil(S)) + 4
    half = N * q + 4 * np.pi * alpha + 20
    eta = np.arange(-half, half + deta, deta)
    K = max(auto_K(S), N) + 3 + int(np.ceil(4 * np.pi * alpha / q)) + Kextra
    nz = required_nz(alpha, aac, S=S, eta_max=half)
    cost_tt = 3e-8 * nz * eta.size
    cost_fl = 2.8e-9 * eta.size * (2 * (2 * K + 1)) ** 3
    use = solver if solver != "auto" else ("tt" if cost_tt < cost_fl else "floquet")
    if use == "tt":
        D0, Dh = solve_laue(alpha, eta, phi=Drive(S=S, alpha_ac=aac).steady_phi(t), n_z=nz)
    else:
        D0, Dh = FloquetLaue(alpha, aac, S, K, t=t).fields(eta)
    return eta, np.abs(Dh) ** 2, np.abs(D0) ** 2, q, N, use


def metrics(alpha, aac, S, nphase=1):
    """nphase > 1: average I over nphase acoustic phases t in [0, 1/aac) (Entry 43b)."""
    Is = [spectrum(alpha, aac, S, t=k / (nphase * aac)) for k in range(nphase)]
    eta, I, I0, q, N, use = Is[0]
    I = np.mean([x[1] for x in Is], axis=0)
    I0 = np.mean([x[2] for x in Is], axis=0)
    ns = np.arange(-N, N + 1)
    a = np.pi * alpha * np.abs(jv(ns, S))
    Iind = sum(laue_line(eta - n * q, an) for n, an in zip(ns, a))
    p = np.array([np.trapezoid(np.where(np.abs(eta - n * q) < q / 2, I, 0), eta) for n in ns])
    pi = np.array([np.trapezoid(np.where(np.abs(eta - n * q) < q / 2, Iind, 0), eta) for n in ns])
    ctr = np.array([np.interp(n * q, eta, I) for n in ns])
    sel = np.abs(jv(ns, S)) > 0.02
    pk = np.flatnonzero((I[1:-1] > I[:-2]) & (I[1:-1] >= I[2:]) & (I[1:-1] > 0.01 * I.max())) + 1
    tot, tot_ind = np.trapezoid(I, eta), np.trapezoid(Iind, eta)
    return dict(alpha=alpha, aac=aac, S=S, R=alpha / aac, solver=use,
                E_spec=float(np.trapezoid(np.abs(I - Iind), eta) / max(tot_ind, 1e-30)) if S > 0 else 0.0,
                E_pop=float(np.abs(p - pi).sum() / max(pi.sum(), 1e-30)),
                E_ctr=float(np.max(np.abs(ctr[sel] - np.sin(a[sel]) ** 2))) if sel.any() else 0.0,
                PR=float(p.sum() ** 2 / max((p ** 2).sum(), 1e-30)),
                PR_ind=float(pi.sum() ** 2 / max((pi ** 2).sum(), 1e-30)),
                tot=float(tot), tot_ind=float(tot_ind), npk=int(pk.size),
                flux=float(np.max(np.abs(I + I0 - 1))),
                p=p.tolist(), p_ind=pi.tolist(), ns=ns.tolist())


NPHASE = int(os.environ.get("E43_NPHASE", "1"))
if NPHASE > 1:
    CACHE = HERE / f"cache_e43_ph{NPHASE}"


def cached(args):
    alpha, aac, S = args
    f = CACHE / f"a{alpha}_q{aac}_S{S}.json"
    if f.exists():
        return json.loads(f.read_text())
    r = metrics(alpha, aac, S, NPHASE)
    f.write_text(json.dumps(r))
    return r


if __name__ == "__main__":
    CACHE.mkdir(exist_ok=True)
    if len(sys.argv) > 1 and sys.argv[1] == "point":
        r = metrics(*map(float, sys.argv[2:5]))
        print({k: v for k, v in r.items() if k not in ("p", "p_ind", "ns")})
        sys.exit()
    grid = [(a, q, s) for s in SS for a in ALPHAS for q in AACS]
    grid.sort(key=lambda x: x[2])
    with Pool(11) as pool:
        for i, r in enumerate(pool.imap_unordered(cached, grid, chunksize=2)):
            if i % 50 == 0:
                print(f"{i}/{len(grid)}", flush=True)
    print("done", flush=True)
