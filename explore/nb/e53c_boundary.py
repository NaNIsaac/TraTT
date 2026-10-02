"""Entry 53c: recompute, with a converged 16-phase average, the Entry 43 points whose 4-phase E_pop lies
within +-30 % of a regime threshold (0.05 / 0.2 / 0.5). Uses e43_phase.metrics unchanged (nphase = 16),
cache in cache_e43_ph16/. Then compares the regime boundaries with the 4-phase ones.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import glob  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import e43_phase  # noqa: E402

C16 = HERE / "cache_e43_ph16"


def run(args):
    a, q, s = args
    f = C16 / f"a{a}_q{q}_S{s}.json"
    if f.exists():
        return json.loads(f.read_text())
    r = e43_phase.metrics(a, q, s, nphase=16)
    f.write_text(json.dumps(r))
    return r


if __name__ == "__main__":
    C16.mkdir(exist_ok=True)
    rows = [json.loads(open(f).read()) for f in glob.glob(str(HERE / "cache_e43_ph4" / "*.json"))]
    sel = [r for r in rows if r["S"] > 0 and any(0.7 * t <= r["E_pop"] <= 1.3 * t for t in (0.05, 0.2, 0.5))]
    print(f"{len(sel)} boundary points of {len(rows)}", flush=True)
    with Pool(11) as pool:
        res = pool.map(run, [(r["alpha"], r["aac"], r["S"]) for r in sel], chunksize=2)
    reg = lambda e: "A" if e < 0.05 else "B" if e < 0.2 else "C" if e < 0.5 else "D"  # noqa: E731
    four = {(r["alpha"], r["aac"], r["S"]): r["E_pop"] for r in sel}
    ch = [(k, four[k], r["E_pop"]) for r in res for k in [(r["alpha"], r["aac"], r["S"])] if reg(four[k]) != reg(r["E_pop"])]
    d = np.array([r["E_pop"] - four[(r["alpha"], r["aac"], r["S"])] for r in res])
    print(f"E_pop(16) - E_pop(4): median {np.median(d):+.4f}, range [{d.min():+.4f}, {d.max():+.4f}]")
    print(f"regime changes: {len(ch)} of {len(res)}")
    for k, e4, e16 in sorted(ch, key=lambda x: x[0][1] and x[0][0] / x[0][1]):
        print(f"  alpha={k[0]} aac={k[1]} S={k[2]} R={k[0]/k[1]:.3f}: {e4:.4f} ({reg(e4)}) -> {e16:.4f} ({reg(e16)})")
