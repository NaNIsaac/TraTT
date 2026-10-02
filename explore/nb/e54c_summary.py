"""Entry 54c: regime boundaries and medians with the exact phase average (cache_e43_exact/) replacing the
16-phase value where available. Same binning as Entry 53 (R bins geomspace(0.003, 5, 30), threshold = first
bin whose median E_pop >= 0.05 / 0.2 / 0.5). Medians printed in scientific notation.
"""
import glob
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(d):
    return {(r["alpha"], r["aac"], r["S"]): r["E_pop"] for r in
            (json.loads(open(f).read()) for f in glob.glob(str(HERE / d / "*.json"))) if r["S"] > 0}


p4, p16, ex = load("cache_e43_ph4"), load("cache_e43_ph16"), load("cache_e43_exact")
merged = {k: p16.get(k, v) for k, v in p4.items()}
final = {k: ex.get(k, v) for k, v in merged.items()}
reg = lambda e: "A" if e < 0.05 else "B" if e < 0.2 else "C" if e < 0.5 else "D"  # noqa: E731
Rb = np.geomspace(0.003, 5, 30)
Rc = np.sqrt(Rb[1:] * Rb[:-1])


def bounds(D, S):
    R, E = np.array([(k[0] / k[1], v) for k, v in D.items() if k[2] == S]).T
    med = np.array([np.median(E[m]) if (m := (R > Rb[i]) & (R <= Rb[i + 1])).any() else np.nan for i in range(29)])
    ok = ~np.isnan(med)
    return " / ".join(f"{Rc[ok][i[0]]:.2f}" if (i := np.flatnonzero(med[ok] >= t)).size else ">5" for t in (0.05, 0.2, 0.5))


if __name__ == "__main__":
    print(f"points: 4ph {len(p4)}, 16ph {len(p16)}, exact {len(ex)}; regime changes merged -> final: "
          f"{sum(reg(merged[k]) != reg(final[k]) for k in ex)}")
    print("boundaries R(0.05 / 0.2 / 0.5): 4ph | merged (16ph) | final (exact)")
    for S in (0.5, 1.0, 2.0, 3.0, 5.0, 8.0):
        print(f"  S={S}: {bounds(p4, S)} | {bounds(merged, S)} | {bounds(final, S)}")
    for lab, sel in (("aac>=2", lambda k: k[1] >= 2), ("aac<2", lambda k: k[1] < 2)):
        for name, D in (("4ph", p4), ("merged", merged), ("final", final)):
            R = np.array([k[0] / k[1] for k in D if sel(k)])
            E = np.array([D[k] for k in D if sel(k)])
            print(f"  {lab:6s} {name:6s}: " + "  ".join(
                f"({lo},{hi}]:{np.median(E[(R > lo) & (R <= hi)]):.1e}" for lo, hi in
                [(0, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 1), (1, 3)]))
    for name, D in (("4ph", p4), ("merged", merged), ("final", final)):
        E = [D[k] for k in D if k[0] <= 0.2]
        print(f"  alpha<=0.2 {name}: median {np.median(E):.1e} max {max(E):.2e}")
