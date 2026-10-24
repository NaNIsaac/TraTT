"""Entry 86 (T2p): does the integrated satellite population lift the beam condition of Entry 84?

Observable: P(t; sigma_eta) = int_channel p(delta) [I(nq + delta, t) - I(nq + delta, 0)] d delta / int_channel p(delta) d delta,
channel |eta - nq| <= q/2, Gaussian beam p of width sigma_eta centred on nq, detector accepting the whole channel.
Reference: the flat-beam limit P_flat = (1/q) int_channel [I(t) - I(0)]. Full TT rocking_cube once per case and geometry
on the channel grid (spacing w_n/40); every sigma_eta is a change of weights.
Pre-registered rule (notebook Entry 86): "integrated observable works with an ordinary beam" if at sigma = 1, 3, 10 w_n
both |dperiod| <= 1 % and |damp| <= 5 % (period = first-maximum time in Laue, half-rise time in Bragg; amplitude =
P / P_flat at that time) in every case and geometry; else the smallest sigma from which both hold.
Usage: python e86_population_divergence.py [run] [fig]
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
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

OUT = HERE / "e86_population_divergence.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e86_population_divergence.png"
CASES = [(0.5, 2.0, 1), (1.5, 2.0, 1), (2.0, 4.0, 2), (1.5, 4.0, 1)]
AAC, OFF = 20.0, 0.0037
T = np.linspace(0.0, 1.0, 61) + OFF
SIG = (0.3, 1.0, 3.0, 10.0)
SCALE_URAD = {"Si440_17.48keV_Laue": 0.12130e-5 / np.sin(np.radians(2 * 21.679)) * 1e6, "Si400_8keV_Bragg": 0.78073e-5 / np.sin(np.radians(2 * 34.802)) * 1e6}


def first_max(t, y):
    i = int(np.argmax(y[: len(y) // 2 + 5]))
    if 0 < i < len(y) - 1:
        a, b, c = y[i - 1], y[i], y[i + 1]
        den = a - 2 * b + c
        dx = 0.5 * (a - c) / den if den != 0 else 0.0
        return t[i] + dx * (t[1] - t[0]), b - 0.25 * (a - c) * dx
    return t[i], y[i]


def half_rise(t, y):
    ym = y[-1]
    j = int(np.argmax(y >= ym / 2))
    if j == 0:
        return t[0], y[0]
    tt = t[j - 1] + (ym / 2 - y[j - 1]) / (y[j] - y[j - 1]) * (t[j] - t[j - 1])
    return tt, ym


def metrics(t, y, geometry):
    return first_max(t, y) if geometry == "laue" else half_rise(t, y)


def run_case(S, alpha, n, geometry, spacing_div=40):
    d = Drive(S=S, alpha_ac=AAC)
    an = alpha * abs(jv(n, S))
    wn = 2 * np.pi * an
    q = d.q
    delta = np.arange(-q / 2, q / 2 + wn / (2 * spacing_div), wn / spacing_div)
    I = rocking_cube(alpha, d, n * q + delta, np.concatenate([[OFF * 0.0 + 1e-6], T]), geometry=geometry, n_z=1200)
    dI = I[1:] - I[0][None, :]  # subtract the unmodulated (t = 0+) rocking curve ray by ray
    return dict(S=S, alpha=alpha, n=n, geometry=geometry, alpha_n=an, w_n=wn, q=q, delta=delta, dI=dI)


def section_run():
    t0 = time.time()
    out = dict(cases=CASES, sigmas=SIG, scale_urad=SCALE_URAD, results=[])
    all_ok = True
    for geometry in ("laue", "bragg"):
        for S, alpha, n in CASES:
            c = run_case(S, alpha, n, geometry)
            delta, dI, wn, q, an = c["delta"], c["dI"], c["w_n"], c["q"], c["alpha_n"]
            dd = delta[1] - delta[0]
            P_flat = dI.sum(axis=1) * dd / q
            t_ref, y_ref = metrics(T, P_flat, geometry)
            rows = []
            for sig in SIG:
                p = np.exp(-0.5 * (delta / (sig * wn)) ** 2)
                p /= p.sum()
                P = (dI * p[None, :]).sum(axis=1)
                tm, ym = metrics(T, P, geometry)
                # amplitude at the reference time: ratio of curves at t_ref (robust to a flat maximum)
                amp = float(np.interp(t_ref, T, P) / np.interp(t_ref, T, P_flat))
                rows.append(dict(sigma=sig, sigma_eta=sig * wn, dperiod=tm / t_ref - 1, damp=amp - 1, P=P.tolist()))
            ok = [r["sigma"] for r in rows if abs(r["dperiod"]) <= 0.01 and abs(r["damp"]) <= 0.05]
            smin = min(ok) if ok else None
            for s_ in (1.0, 3.0, 10.0):
                r = [x for x in rows if x["sigma"] == s_][0]
                if abs(r["dperiod"]) > 0.01 or abs(r["damp"]) > 0.05:
                    all_ok = False
            res = dict(geometry=geometry, S=S, alpha=alpha, n=n, alpha_n=an, w_n=wn, w_n_hw=abs(jv(n, S)), q_hw=AAC / alpha, n_rays=int(delta.size),
                       P_flat=P_flat.tolist(), t_ref=t_ref, rows=rows, sigma_min_ok=smin, P_flat_end_over_centre=float(P_flat[-1]))
            if geometry == "laue":
                law = np.array([quad(j0, 0, 2 * np.pi * an * tt)[0] for tt in T])
                res["law_shape_err"] = float(np.max(np.abs(P_flat / P_flat[-1] - law / law[-1])))
            else:
                law = np.tanh(np.pi * an * T) ** 2
                res["law_shape_err"] = float(np.max(np.abs(P_flat / P_flat[-1] - law / law[-1])))
            out["results"].append(res)
            print(f"{geometry:5s} S={S:g} a={alpha:g} n={n} ({delta.size} rays, w_n = {abs(jv(n, S)):.3f} hw, channel q = {AAC / alpha:.1f} hw): "
                  + " | ".join(f"s={r['sigma']:g}: dT {r['dperiod'] * 100:+.2f} %, dA {r['damp'] * 100:+.1f} %" for r in rows)
                  + f" | smallest ok sigma {smin} | P_flat shape vs law max dev {res['law_shape_err']:.3f} | P_flat(1) = {P_flat[-1]:.4f}", flush=True)
    # convergence: spacing halved, one case
    c1 = run_case(1.5, 2.0, 1, "laue", spacing_div=80)
    p = np.exp(-0.5 * (c1["delta"] / (3 * c1["w_n"])) ** 2)
    p /= p.sum()
    P80 = (c1["dI"] * p[None, :]).sum(axis=1)
    r = [x for x in out["results"] if x["geometry"] == "laue" and x["S"] == 1.5 and x["alpha"] == 2.0][0]
    P40 = np.array([x for x in r["rows"] if x["sigma"] == 3.0][0]["P"])
    conv = float(np.max(np.abs(P80 / P40 - 1)))
    out["conv_spacing"] = conv
    print(f"convergence (laue S=1.5 a=2, sigma = 3 w_n, spacing w_n/40 -> w_n/80): max |dP/P| = {conv * 100:.3f} %")
    out["verdict"] = "integrated observable works with an ordinary beam (sigma = 1, 3, 10 w_n all within 1 % / 5 %)" if all_ok else "condition table (smallest ok sigma per case)"
    print("verdict:", out["verdict"], f"({time.time() - t0:.0f} s)")
    json.dump(out, open(OUT, "w"), indent=1)


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, geometry in zip(axs, ("laue", "bragg"), strict=True):
        r = [x for x in d["results"] if x["geometry"] == geometry and x["S"] == 1.5 and x["alpha"] == 2.0][0]
        an = r["alpha_n"]
        Pf = np.array(r["P_flat"])
        ax.plot(an * T, Pf / Pf[-1], "k-", lw=2, label="flat beam (channel integral)")
        for row, c in zip(r["rows"], ("C0", "C1", "C3", "0.6"), strict=True):
            P = np.array(row["P"])
            ax.plot(an * T, P / Pf[-1], "--", color=c, lw=1.3, label=f"$\\sigma_\\eta$ = {row['sigma']:g} $w_n$")
        ax.set(xlabel="$\\alpha|J_n(S)|\\,\\tilde t$", ylabel="integrated population / flat-beam value at $\\tilde t$ = 1", title=f"{geometry.capitalize()}, S = 1.5, $\\alpha$ = 2, n = 1")
        ax.grid(color="#e4e3df", lw=0.6)
        ax.legend(fontsize=8, frameon=False)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
