"""Entry 53b: does the converged phase average change the Entry 43 regime metrics?

For representative points of the phase-averaged cache (explore/nb/cache_e43_ph4, 4 phases), recompute
the order-window population error E_pop with the exact phase average <I(eta)> = sum_m |b_m(eta)|^2
(Entry 53), using the same eta grid, windows and independent-channel reference as
e43_phase.metrics, and compare. Also reports the 16-phase TT average as an independent check.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.special import jv  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e43_phase import laue_line, spectrum  # noqa: E402
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402

# (alpha, alpha_ac, S): boundaries R ~ 0.3 / 0.75 / 2, small and large S, integer and non-integer alpha_ac
POINTS = [(1.0, 3.0, 1.0), (1.0, 3.0, 5.0), (0.5, 1.5, 1.0), (0.5, 1.5, 5.0), (0.2, 0.75, 1.0),
          (2.0, 3.0, 1.0), (5.0, 8.0, 1.0), (5.0, 8.0, 5.0), (1.0, 1.5, 1.0), (1.0, 1.5, 5.0), (0.5, 0.75, 1.0),
          (5.0, 2.0, 1.0), (5.0, 2.0, 5.0), (3.0, 1.5, 1.0), (3.0, 1.5, 5.0), (1.0, 0.5, 1.0)]


def epop(eta, I, alpha, aac, S):
    q = 2 * np.pi * aac
    N = int(np.ceil(S)) + 4
    ns = np.arange(-N, N + 1)
    Iind = sum(laue_line(eta - n * q, np.pi * alpha * abs(jv(n, S))) for n in ns)
    p = np.array([np.trapezoid(np.where(np.abs(eta - n * q) < q / 2, I, 0), eta) for n in ns])
    pi = np.array([np.trapezoid(np.where(np.abs(eta - n * q) < q / 2, Iind, 0), eta) for n in ns])
    return float(np.abs(p - pi).sum() / pi.sum())


def regime(e):
    return "A" if e < 0.05 else "B" if e < 0.2 else "C" if e < 0.5 else "D"


if __name__ == "__main__":
    out = []
    print(" alpha  aac   S    R    | E_pop cache(4ph) | 4ph recomputed | 16ph TT | exact avg | regime 4ph -> exact")
    for alpha, aac, S in POINTS:
        q = 2 * np.pi * aac
        cache = json.loads((HERE / "cache_e43_ph4" / f"a{alpha}_q{aac}_S{S}.json").read_text())
        eta = spectrum(alpha, aac, S, solver="tt")[0]
        I4 = np.mean([spectrum(alpha, aac, S, t=k / (4 * aac), solver="tt")[1] for k in range(4)], axis=0)
        I16 = np.mean([spectrum(alpha, aac, S, t=k / (16 * aac), solver="tt")[1] for k in range(16)], axis=0)
        N = int(np.ceil(S)) + 4
        K = max(auto_K(S), N) + 3 + int(np.ceil(4 * np.pi * alpha / q))
        x = FloquetLaue(alpha, aac, S, K).propagate(eta)
        Iex = np.sum(np.abs(x[:, x.shape[1] // 2:]) ** 2, axis=1)
        e4, e16, ex = epop(eta, I4, alpha, aac, S), epop(eta, I16, alpha, aac, S), epop(eta, Iex, alpha, aac, S)
        out.append(dict(alpha=alpha, aac=aac, S=S, R=alpha / aac, cache=cache["E_pop"], e4=e4, e16=e16, exact=ex,
                        maxdiff_I16_exact=float(np.abs(I16 - Iex).max())))
        print(f" {alpha:4.1f} {aac:5.2f} {S:3.1f} {alpha/aac:5.2f} | {cache['E_pop']:.4f} | {e4:.4f} | {e16:.4f} | "
              f"{ex:.4f} | {regime(e4)} -> {regime(ex)}  (max|I16 - Iexact| {out[-1]['maxdiff_I16_exact']:.1e})", flush=True)
    json.dump(out, open(HERE / "e53b_map_check.json", "w"), indent=1)
