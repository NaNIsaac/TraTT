"""Entry 79 (T2e): does the coupling ratio R = alpha / alpha_ac break the running-mean degeneracy class?

Entry 76: at R = 0.0125 two histories with the same running mean of J_n(S A) give the same order-n satellite
to 3e-6 for S <= 6; one exploratory row at R = 0.1 gave 2.9e-4. Here R is scanned at fixed alpha_ac = 10
(alpha = 0.125 ... 3), constructions (a) S = 3 and (a2) S = 6 of Entry 76.

Pre-registered rule (notebook Entry 79): "R-boundary exists" if at some R <= 0.3 the degenerate-order TT
difference in the window exceeds 5e-4 (or 10x the measured convergence error), the smallest such R being
the boundary, reported with the law's own error there; "R-independent" otherwise.
Usage: python e79_degeneracy_R.py [run] [fig]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e76_degeneracy_S import one_case  # noqa: E402

OUT = HERE / "e79_degeneracy_R.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e79_degeneracy_R.png"
AAC = 10.0
ALPHAS = (0.125, 0.5, 1.0, 2.0, 3.0)


def section_run():
    out = dict(aac=AAC, alphas=ALPHAS, cases=[])
    t0 = time.time()
    for alpha in ALPHAS:
        for S, key in ((3.0, "a"), (6.0, "a2")):
            r = one_case(S, key, alpha=alpha, aac=AAC)
            r["R"] = alpha / AAC
            out["cases"].append(r)
    bnd = {}
    for key in ("a", "a2"):
        rows = [r for r in out["cases"] if r["key"] == key]
        R = np.array([r["R"] for r in rows])
        d = np.array([r["n1"]["dtt_win"] for r in rows])
        err = np.array([r["n1"]["law_err"] for r in rows])
        hit = [r for r in rows if r["broken"]]
        slope = np.polyfit(np.log(R), np.log(d), 1)[0]
        bnd[key] = dict(
            R=R.tolist(),
            d=d.tolist(),
            law_err=err.tolist(),
            below_law_err=bool(np.all(d < err)),
            slope=float(slope),
            boundary=(hit[0]["R"] if hit else None),
            boundary_law_err=(hit[0]["n1"]["law_err"] if hit else None),
        )
        print(
            f"construction {key}: residual vs R = "
            + " ".join(f"{a:g}:{b:.1e}" for a, b in zip(R, d, strict=True))
            + f"; law error {' '.join(f'{v:.1e}' for v in err)}; power law R^{slope:.2f}; "
            + (f"R-boundary at R = {hit[0]['R']:g} (law error there {hit[0]['n1']['law_err']:.1e})" if hit else "no boundary for R <= 0.3")
        )
    out["bounds"] = bnd
    anyb = [b["boundary"] for b in bnd.values() if b["boundary"] is not None]
    out["verdict"] = f"R-boundary exists, smallest R = {min(anyb):g}" if anyb else "R-independent for R <= 0.3"
    print("verdict:", out["verdict"], f"({time.time() - t0:.0f} s)")
    json.dump(out, open(OUT, "w"), indent=1)


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    fig, ax = plt.subplots(figsize=(5.4, 4), constrained_layout=True)
    for key, mk, lab in (("a", "o", "S = 3, J$_1$ rising"), ("a2", "s", "S = 6, J$_1$ descending")):
        b = d["bounds"][key]
        (ln,) = ax.plot(b["R"], b["d"], mk + "-", label=f"{lab}: TT $|\\Delta I_1|$, same $\\Theta_1$")
        ax.plot(b["R"], b["law_err"], mk + ":", color=ln.get_color(), mfc="none", label=f"{lab}: law error")
    ax.axhline(5e-4, color="C3", lw=1, ls="--")
    ax.text(0.013, 6e-4, "threshold 5e-4", color="C3", fontsize=8)
    ax.set(xscale="log", yscale="log", xlabel="coupling ratio $R = \\alpha/\\alpha_{ac}$", ylabel="max over the window", title="degeneracy residual against R ($\\alpha_{ac}$ = 10)")
    ax.legend(fontsize=7)
    ax.grid(color="#e4e3df", lw=0.6)
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
