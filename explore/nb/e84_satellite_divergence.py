"""Entry 84 (T2n): beam divergence and the formation law at the satellite centre (the poster's observable, P2).

Four poster cases (S, alpha, n), alpha_ac = 20, step front, Laue and Bragg; full TT rocking_cube on eta = n q + delta,
Gaussian beam weights of width sigma_eta = (0.1, 0.3, 1) x w_n with w_n = 2 pi alpha |J_n(S)| (the satellite's own
Darwin-like half-width). Metrics on the first maximum of I_n(t): its time and height, beam-averaged vs single ray.
Pre-registered rule (notebook Entry 84): "P2 measurable with an ordinary beam" if at sigma = 1 w_n |dperiod| <= 1 % and
|damp| <= 5 % in every case and geometry; else the largest sigma meeting both. Check: the beam-averaged
independent-satellite closed form (lossless two-beam slab of thickness t~ with coupling alpha_n at detuning delta, via
rocking(alpha_n t~, delta t~)) within 2 % of the TT beam average. Reported: sigma = 3 w_n vs the integrated-population
law int_0^{2 pi alpha_n t~} J_0.
Usage: python e84_satellite_divergence.py [run] [fig]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import j0, jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.core import rocking  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

OUT = HERE / "e84_satellite_divergence.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e84_satellite_divergence.png"
CASES = [(0.5, 2.0, 1), (1.5, 2.0, 1), (2.0, 4.0, 2), (1.5, 4.0, 1)]
AAC, OFF = 20.0, 0.0037
T = np.linspace(0.0, 1.0, 61) + OFF
SIG = (0.1, 0.3, 1.0, 3.0)
# Darwin half-widths in angle (|chi_rh| / sin 2 theta_B), X0h constants of e20_units_x0h.py
SCALE = {"Si440_17.48keV_Laue": 0.12130e-5 / np.sin(np.radians(2 * 21.679)), "Si400_8keV_Bragg": 0.78073e-5 / np.sin(np.radians(2 * 34.802))}


def first_max(t, y):
    i = int(np.argmax(y[: len(y) // 2 + 5]))  # first maximum lies in the first half period for alpha_n t <= 0.5
    if 0 < i < len(y) - 1:
        a, b, c = y[i - 1], y[i], y[i + 1]
        den = a - 2 * b + c
        dx = 0.5 * (a - c) / den if den != 0 else 0.0
        return t[i] + dx * (t[1] - t[0]), b - 0.25 * (a - c) * dx
    return t[i], y[i]


def closed_form(alpha_n, delta, t, geometry):
    """Independent-satellite two-beam slab: thickness t~ scales alpha -> alpha_n t~, eta -> delta t~."""
    out = np.empty((t.size, delta.size))
    for i, tt in enumerate(t):
        out[i] = rocking(alpha_n * tt, delta * tt, phi=None, geometry=geometry, n_z=512)
    return out


def section_run():
    t0 = time.time()
    out = dict(cases=CASES, sigmas=SIG, scale_urad={k: v * 1e6 for k, v in SCALE.items()}, results=[])
    verdict_ok = True
    for geometry in ("laue", "bragg"):
        for S, alpha, n in CASES:
            d = Drive(S=S, alpha_ac=AAC)
            an = alpha * abs(jv(n, S))
            wn = 2 * np.pi * an
            tmax_s = None
            rows = []
            for sig in (0.0,) + SIG:
                se = sig * wn
                if se == 0:
                    delta, p = np.array([0.0]), np.array([1.0])
                else:
                    delta = np.arange(-5 * se, 5 * se + wn / 80, wn / 40)
                    p = np.exp(-0.5 * (delta / se) ** 2)
                    p /= p.sum()
                I = rocking_cube(alpha, d, n * d.q + delta, T, geometry=geometry, n_z=1200)
                Ic = (I * p[None, :]).sum(axis=1)
                tm, Im = first_max(T, Ic)
                C = closed_form(an, delta, T, geometry)
                Cc = (C * p[None, :]).sum(axis=1)
                cf_err = float(np.max(np.abs(Ic - Cc)))
                if sig == 0:
                    tmax_s, Imax_s, I_single = tm, Im, Ic
                row = dict(sigma=sig, sigma_eta=se, tmax=tm, Imax=Im, dperiod=tm / tmax_s - 1, damp=Im / Imax_s - 1, cf_abs_err=cf_err, I=Ic.tolist())
                if sig == 3.0 and geometry == "laue":
                    pop = np.array([quad(j0, 0, 2 * np.pi * an * tt)[0] for tt in T])
                    k = Ic[-1] / pop[-1]
                    row["pop_law_err"] = float(np.max(np.abs(Ic - k * pop)) / Ic.max())
                rows.append(row)
                print(f"{geometry:5s} S={S:g} a={alpha:g} n={n}: sigma={sig:3g} w_n ({se:6.2f} eta): t_max {tm:.4f} ({row['dperiod'] * 100:+.2f} %), I_max {Im:.4f} ({row['damp'] * 100:+.1f} %), TT-vs-CF max|dI| {cf_err:.1e}"
                      + (f", population law err {row['pop_law_err']:.3f}" if "pop_law_err" in row else ""), flush=True)
            r1 = [r for r in rows if r["sigma"] == 1.0][0]
            if abs(r1["dperiod"]) > 0.01 or abs(r1["damp"]) > 0.05:
                verdict_ok = False
            ok = [r["sigma"] for r in rows if r["sigma"] > 0 and abs(r["dperiod"]) <= 0.01 and abs(r["damp"]) <= 0.05]
            smax = max(ok) if ok else 0.0
            out["results"].append(dict(geometry=geometry, S=S, alpha=alpha, n=n, alpha_n=an, w_n=wn, w_n_hw=abs(jv(n, S)), sigma_max_wn=smax, rows=rows))
            print(f"   largest sigma meeting both (scan): {smax:g} w_n = {smax * abs(jv(n, S)):.3f} main half-widths; w_n = {abs(jv(n, S)):.3f} half-widths = "
                  + ", ".join(f"{abs(jv(n, S)) * v * 1e6:.2f} urad ({k})" for k, v in SCALE.items()))
    out["verdict"] = "P2 measurable with an ordinary beam" if verdict_ok else "condition table applies (sigma_eta < w_n)"
    print("verdict:", out["verdict"], f"({time.time() - t0:.0f} s)")
    json.dump(out, open(OUT, "w"), indent=1)


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True, sharey=True)
    for ax, geometry in zip(axs, ("laue", "bragg"), strict=True):
        r = [x for x in d["results"] if x["geometry"] == geometry and x["S"] == 1.5 and x["alpha"] == 2.0][0]
        an = r["alpha_n"]
        for row, c in zip(r["rows"], ("k", "C0", "C1", "C3", "0.6"), strict=True):
            ax.plot(an * T, row["I"], "-", color=c, lw=1.4, label=f"$\\sigma_\\eta$ = {row['sigma']:g} $w_n$")
        ax.set(xlabel="$\\alpha|J_n(S)|\\,\\tilde t$", title=f"{geometry.capitalize()}: S = 1.5, $\\alpha$ = 2, n = 1")
        ax.grid(color="#e4e3df", lw=0.6)
        ax.legend(fontsize=8, frameon=False)
    axs[0].set_ylabel("$I_n$ (beam average)")
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
