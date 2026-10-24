"""Entry 87 (T2q): the Laue integrated satellite population at sigma_eta = 3 w_n, corrected metrics.

Definitions (notebook Entry 87, fixed before the run):
    C_n = {eta : |eta - n q| <= q/2},  q = 2 pi alpha_ac;  dI(eta, t) = I(eta, t) - I(eta, 0+)   (full TT rocking_cube)
    P_sigma(t) = int_C p_sigma(delta) dI(nq + delta, t) d delta / int_C p_sigma d delta,  p_sigma Gaussian, sigma_eta = 3 w_n
    P_0(t)     = (1/q) int_C dI d delta                                                  (flat beam; the law)
    S^_sigma(t) = P_sigma(t) / P_sigma(t_end),  S^_0(t) = P_0(t) / P_0(t_end)              (t_end = 1 t_tr step, 2 t_tr ramp)
    t*         = time of the first maximum of S^ (parabolic interpolation)
    dperiod    = t*_sigma / t*_0 - 1;   A = S^_sigma(t*_0) / S^_0(t*_0) - 1
Rule: |dperiod| <= 1 % and |A| <= 5 % in all five runs (four step cases + one ramp, rise = 1 t_tr) -> "measurable in the
high-resolution beamline range"; else the condition table (sigma = 3, 5, 10 w_n scanned, smallest passing).
Usage: python e87_population_laue.py [run] [fig]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

OUT = HERE / "e87_population_laue.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e87_population_laue.png"
AAC, OFF = 20.0, 0.0037
RUNS = [(0.5, 2.0, 1, 0.0), (1.5, 2.0, 1, 0.0), (2.0, 4.0, 2, 0.0), (1.5, 4.0, 1, 0.0), (1.5, 2.0, 1, 1.0)]  # (S, alpha, n, rise)
SIG = (3.0, 5.0, 10.0)
URAD_PER_HW = 0.12130e-5 / np.sin(np.radians(2 * 21.679)) * 1e6  # Si 440 / 17.48 keV Darwin half-width


def first_max(t, y):
    i = int(np.argmax(y[: int(len(y) * 0.6)]))
    if 0 < i < len(y) - 1:
        a, b, c = y[i - 1], y[i], y[i + 1]
        den = a - 2 * b + c
        dx = 0.5 * (a - c) / den if den != 0 else 0.0
        return t[i] + dx * (t[1] - t[0])
    return t[i]


def run_case(S, alpha, n, rise, spacing_div=40):
    d = Drive(S=S, alpha_ac=AAC, rise=rise)
    an = alpha * abs(jv(n, S))
    wn = 2 * np.pi * an
    q = d.q
    t_end = 1.0 if rise == 0 else 2.0
    T = np.linspace(0.0, t_end, 61 if rise == 0 else 121) + OFF
    delta = np.arange(-q / 2, q / 2 + wn / (2 * spacing_div), wn / spacing_div)
    I = rocking_cube(alpha, d, n * q + delta, np.concatenate([[1e-6], T]), geometry="laue", n_z=1200)
    dI = I[1:] - I[0][None, :]
    return dict(S=S, alpha=alpha, n=n, rise=rise, alpha_n=an, w_n=wn, q=q, T=T, delta=delta, dI=dI)


def shapes(c, sig):
    delta, dI, wn, q = c["delta"], c["dI"], c["w_n"], c["q"]
    dd = delta[1] - delta[0]
    P0 = dI.sum(axis=1) * dd / q
    if sig == 0:
        P = P0
    else:
        p = np.exp(-0.5 * (delta / (sig * wn)) ** 2)
        P = (dI * p[None, :]).sum(axis=1) / p.sum()
    return P / P[-1], P0 / P0[-1], P, P0


def section_run():
    t0 = time.time()
    out = dict(runs=RUNS, sigmas=SIG, urad_per_hw=URAD_PER_HW, results=[])
    all_ok = True
    for S, alpha, n, rise in RUNS:
        c = run_case(S, alpha, n, rise)
        T = c["T"]
        rows = []
        for sig in SIG:
            Ss, S0, P, P0 = shapes(c, sig)
            t0s, t00 = first_max(T, Ss), first_max(T, S0)
            dperiod = t0s / t00 - 1
            A = float(np.interp(t00, T, Ss) / np.interp(t00, T, S0) - 1)
            rows.append(dict(sigma=sig, sigma_eta=sig * c["w_n"], sigma_hw=sig * abs(jv(n, S)), sigma_urad=sig * abs(jv(n, S)) * URAD_PER_HW,
                             t_star_0=t00, t_star_sigma=t0s, dperiod=dperiod, A=A, ok=bool(abs(dperiod) <= 0.01 and abs(A) <= 0.05), shape=Ss.tolist()))
            if sig == 3.0 and not rows[-1]["ok"]:
                all_ok = False
        ok = [r["sigma"] for r in rows if r["ok"]]
        res = dict(S=S, alpha=alpha, n=n, rise=rise, alpha_n=c["alpha_n"], w_n_hw=abs(jv(n, S)), q_hw=AAC / alpha, q_urad=AAC / alpha * URAD_PER_HW,
                   n_rays=int(c["delta"].size), T=T.tolist(), P0=shapes(c, 0)[3].tolist(), law_shape=shapes(c, 0)[1].tolist(), rows=rows, smallest_ok=(min(ok) if ok else None))
        out["results"].append(res)
        print(f"S={S:g} a={alpha:g} n={n} rise={rise:g} ({c['delta'].size} rays, w_n = {abs(jv(n, S)):.3f} hw, q = {AAC / alpha:.0f} hw = {AAC / alpha * URAD_PER_HW:.0f} urad): "
              + " | ".join(f"s={r['sigma']:g} w_n ({r['sigma_urad']:.2f} urad): dperiod {r['dperiod'] * 100:+.2f} %, A {r['A'] * 100:+.2f} % {'ok' if r['ok'] else 'FAIL'}" for r in rows)
              + f" | smallest ok {res['smallest_ok']}", flush=True)
    c = run_case(1.5, 2.0, 1, 0.0, spacing_div=80)
    Ss80 = shapes(c, 3.0)[0]
    r = [x for x in out["results"] if x["S"] == 1.5 and x["alpha"] == 2.0 and x["rise"] == 0.0][0]
    Ss40 = np.array([x for x in r["rows"] if x["sigma"] == 3.0][0]["shape"])
    conv = float(np.max(np.abs(Ss80 - Ss40)))
    out["conv_spacing"] = conv
    print(f"convergence (S=1.5 a=2 step, sigma = 3 w_n, w_n/40 -> w_n/80): max |dS^| = {conv:.2e}")
    out["verdict"] = ("P2 integrated observable measurable in the high-resolution beamline range (sigma_eta = 3 w_n, all five runs within 1 % / 5 %)"
                      if all_ok else "condition table (smallest passing sigma per run)")
    print("verdict:", out["verdict"], f"({time.time() - t0:.0f} s)")
    json.dump(out, open(OUT, "w"), indent=1)


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, (S, alpha, n, rise) in zip(axs, ((1.5, 2.0, 1, 0.0), (1.5, 2.0, 1, 1.0)), strict=True):
        r = [x for x in d["results"] if x["S"] == S and x["alpha"] == alpha and x["rise"] == rise][0]
        T = np.array(r["T"])
        an = r["alpha_n"]
        ax.plot(an * T, r["law_shape"], "k-", lw=2, label="flat beam $\\hat S_0$ (the law: TT channel integral)")
        for row, c in zip(r["rows"], ("C0", "C1", "C3"), strict=True):
            ax.plot(an * T, row["shape"], "--", color=c, lw=1.3, label=f"$\\sigma_\\eta$ = {row['sigma']:g} $w_n$ = {row['sigma_urad']:.1f} µrad")
        ax.plot(an * T, np.sin(np.pi * an * T) ** 2 * max(r["law_shape"]) if rise == 0 else np.full(T.size, np.nan), ":", color="0.5", lw=1, label="centre law sin^2, scaled to the same maximum" if rise == 0 else None)
        ax.set(xlabel="$\\alpha|J_n(S)|\\,\\tilde t$", ylabel="normalised integrated population $\\hat S$", title=f"Laue, S = {S:g}, $\\alpha$ = {alpha:g}, n = {n}, {'step' if rise == 0 else 'ramp 1 t_tr'}")
        ax.grid(color="#e4e3df", lw=0.6)
        ax.legend(fontsize=7, frameon=False)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
