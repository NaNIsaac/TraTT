"""Entry 54a: exact phase-averaged map points, <I(eta)> = sum_m |b_m(eta)|^2 (Entry 53 identity), with the
same eta grid, windows and independent reference as e43_phase.metrics. Targets: all S = 5 and 8 points,
and points whose merged (16-phase where available) E_pop lies within +-30 % of 0.05 / 0.2 / 0.5.
Cache: cache_e43_exact/. Compares with the 16-phase cache (cache_e43_ph16/).
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import glob  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.special import jv  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e43_phase import laue_line, spectrum  # noqa: E402
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402

CEX = HERE / "cache_e43_exact"


def exact_point(args):
    alpha, aac, S = args
    f = CEX / f"a{alpha}_q{aac}_S{S}.json"
    if f.exists():
        return json.loads(f.read_text())
    q = 2 * np.pi * aac
    N = int(np.ceil(S)) + 4
    eta = spectrum(alpha, aac, S, solver="tt")[0]          # identical eta grid to e43
    K = max(auto_K(S), N) + 3 + int(np.ceil(4 * np.pi * alpha / q))
    out = {}
    for KK in (K, K + 4):                                   # truncation check
        x = FloquetLaue(alpha, aac, S, KK).propagate(eta)
        I = np.sum(np.abs(x[:, x.shape[1] // 2:]) ** 2, axis=1)
        ns = np.arange(-N, N + 1)
        Iind = sum(laue_line(eta - n * q, np.pi * alpha * abs(jv(n, S))) for n in ns)
        p = np.array([np.trapezoid(np.where(np.abs(eta - n * q) < q / 2, I, 0), eta) for n in ns])
        pi = np.array([np.trapezoid(np.where(np.abs(eta - n * q) < q / 2, Iind, 0), eta) for n in ns])
        out[KK] = float(np.abs(p - pi).sum() / pi.sum())
    r = dict(alpha=alpha, aac=aac, S=S, E_pop=out[K], E_pop_K4=out[K + 4], K=K)
    f.write_text(json.dumps(r))
    return r


if __name__ == "__main__":
    CEX.mkdir(exist_ok=True)

    def load(d):
        return {(r["alpha"], r["aac"], r["S"]): r for r in
                (json.loads(open(f).read()) for f in glob.glob(str(HERE / d / "*.json")))}

    p4, p16 = load("cache_e43_ph4"), load("cache_e43_ph16")
    merged = {k: (p16[k] if k in p16 else v)["E_pop"] for k, v in p4.items() if k[2] > 0}
    todo = sorted({k for k, e in merged.items()
                   if k[2] in (5.0, 8.0) or any(0.7 * t <= e <= 1.3 * t for t in (0.05, 0.2, 0.5))})
    print(f"{len(todo)} points", flush=True)
    with Pool(11) as pool:
        res = pool.map(exact_point, todo, chunksize=1)
    json.dump(res, open(HERE / "e54_exactmap.json", "w"))
    print("done", flush=True)
