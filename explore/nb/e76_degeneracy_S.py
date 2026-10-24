"""Entry 76 (T2b): does the running-mean degeneracy class of Entry 72 survive larger S?

Entry 72: order n of the transient reads only Theta_n(t) = int_0^1 J_n(S A(t - zeta)) d zeta (Laue: I_n =
sin^2(pi alpha Theta_n)), so two histories with the same running mean of J_n(S A) are degenerate at order n.
Confirmed at S = 1.5 (TT difference 9e-7). Here the same construction is repeated for S = 1.5-6, at both
orders, with the wobble placed on a monotone branch of J_n(S s) so that the inversion exists:
    history 1: A_1(s) = clip(s, 0, 1);   history 2: J_n(S A_2(s)) = J_n(S A_1(s)) + delta(s),
    delta = eps sin(2 pi (s - a)/(b - a)) on [a, b], zero outside; [a, b] inside one monotone branch.
For t in [b, a + 1] the window [t - 1, t] holds the whole wobble, so Theta_n is identical there.

Pre-registered rule (notebook Entry 76): Laue, alpha = 0.5, alpha_ac = 40 (R = 0.0125), 16 carrier phases,
eta = q and 2q, 151 times on [0, 1.5]. "Broken at order n at this S" if max_window |I_n^TT(2) - I_n^TT(1)|
> 5e-4 (ten times the Entry 72 convergence error; if the re-measured convergence error at some S exceeds
5e-5 the threshold there is ten times that and is reported). "Class holds" if no S <= 6 exceeds it.
Side checks: non-degenerate order TT vs law difference to 1e-3; |I_TT - I_law| per S.
Exploratory, not under the rule: construction (a) at S = 3 with R = 0.1 (alpha = 1, alpha_ac = 10).
Usage: python e76_degeneracy_S.py [run] [fig]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.core import required_nz, rocking  # noqa: E402

OUT = HERE / "e76_degeneracy_S.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e76_degeneracy_S.png"
TIMES = np.linspace(0.0, 1.5, 151)
S_SCAN = (1.5, 2.0, 3.0, 4.0, 5.0, 6.0)
THRESH = 5e-4
# monotone branches of J_n: (order, x_lo, x_hi, eps, label)
BRANCHES = {
    "a": (1, 0.0, 1.8412, 0.1, "J1 rising"),
    "a2": (1, 1.8412, 5.3314, 0.1, "J1 descending"),
    "b": (2, 0.0, 3.0542, 0.05, "J2 rising"),
}


def A1(s):
    return np.clip(s, 0.0, 1.0)


def choose_interval(S, n, x_lo, x_hi, eps):
    """[a, b] in s such that J_n(S s) +- eps stays inside the branch range, 0.02 <= a < b <= 0.95."""
    x = np.linspace(x_lo, x_hi, 20001)
    j = jv(n, x)
    lo, hi = j[0], j[-1]
    ok = (np.abs(j - lo) >= eps) & (np.abs(hi - j) >= eps)
    if not ok.any():
        return None
    xa, xb = x[ok][0], x[ok][-1]
    a, b = max(xa / S, 0.02), min(xb / S, 0.95)
    return (a, b) if b - a > 0.1 else None


def make_history2(S, n, x_lo, x_hi, eps, a, b, ns=30001):
    s = np.linspace(0.0, 1.5, ns)
    delta = np.where((s >= a) & (s <= b), eps * np.sin(2 * np.pi * (s - a) / (b - a)), 0.0)
    y = jv(n, S * A1(s)) + delta
    A2 = A1(s).copy()
    for i in np.nonzero(delta != 0)[0]:
        A2[i] = brentq(lambda x, yy=y[i]: jv(n, x) - yy, x_lo, x_hi, xtol=1e-13) / S
    return lambda q: np.interp(q, s, A2, left=0.0, right=A2[-1])


def theta(env, S, n, t, nz=8001):
    z = np.linspace(0.0, 1.0, nz)
    return np.array([np.trapezoid(jv(n, S * env(x - z)), z) for x in np.atleast_1d(t)])


def law(env, S, n, t, alpha):
    return np.sin(np.pi * alpha * theta(env, S, n, t)) ** 2


def tt_trace(env, S, alpha, aac, times, phases=16, nz=None):
    """Phase-averaged I at eta = q and 2q; returns (n_t, 2)."""
    Q = 2 * np.pi * aac
    nz = nz or required_nz(alpha, aac, S=S, eta_max=2 * Q)
    tt = np.asarray(times, float)[:, None, None]
    pp = (np.arange(phases) * 2 * np.pi / phases)[None, :, None]

    def phi(z):
        return S * env(tt - z) * np.sin(Q * (z - tt) + pp)

    inten = rocking(alpha, np.array([Q, 2 * Q]), phi=phi, n_z=nz)
    return inten.mean(axis=1)


def one_case(S, key, alpha=0.5, aac=40.0, verbose=True):
    n, x_lo, x_hi, eps, label = BRANCHES[key]
    iv = choose_interval(S, n, x_lo, x_hi, eps)
    if iv is None:
        return None
    a, b = iv
    env2 = make_history2(S, n, x_lo, x_hi, eps, a, b)
    win = (TIMES >= b) & (TIMES <= min(a + 1.0, TIMES[-1]))
    t0 = time.time()
    res = dict(S=S, key=key, n=n, label=label, a=a, b=b, eps=eps, alpha=alpha, aac=aac)
    ss = np.linspace(0, 1.5, 3001)
    res["maxdA"] = float(np.max(np.abs(env2(ss) - A1(ss))))
    L = {m: (law(A1, S, m, TIMES, alpha), law(env2, S, m, TIMES, alpha)) for m in (1, 2)}
    T1, T2 = tt_trace(A1, S, alpha, aac, TIMES), tt_trace(env2, S, alpha, aac, TIMES)
    for m in (1, 2):
        d_law = L[m][1] - L[m][0]
        d_tt = T2[:, m - 1] - T1[:, m - 1]
        res[f"n{m}"] = dict(
            dlaw_win=float(np.max(np.abs(d_law[win]))),
            dtt_win=float(np.max(np.abs(d_tt[win]))),
            dlaw_out=float(np.max(np.abs(d_law[~win]))),
            dtt_out=float(np.max(np.abs(d_tt[~win]))),
            d_tt_vs_law=float(np.max(np.abs(d_tt - d_law))),
            law_err=float(max(np.max(np.abs(T1[:, m - 1] - L[m][0])), np.max(np.abs(T2[:, m - 1] - L[m][1])))),
            I_max=float(max(T1[:, m - 1].max(), T2[:, m - 1].max())),
        )
    sub = np.array([0.5, 0.95, 1.4])
    nz = required_nz(alpha, aac, S=S, eta_max=4 * np.pi * aac)
    c1 = tt_trace(env2, S, alpha, aac, sub, phases=16, nz=nz)
    c2 = tt_trace(env2, S, alpha, aac, sub, phases=32, nz=2 * nz)
    res["conv"] = [float(v) for v in np.max(np.abs(c1 - c2), axis=0)]
    res["nz"] = int(nz)
    res["thresh"] = float(max(THRESH, 10 * res["conv"][n - 1]))
    res["broken"] = bool(res[f"n{n}"]["dtt_win"] > res["thresh"])
    res["time_s"] = time.time() - t0
    if verbose:
        o = 3 - n
        print(
            f"S={S:3g} {key:>2} ({label:13s}) [a,b]=[{a:.3f},{b:.3f}] max|dA|={res['maxdA']:.3f} n_z={nz}: "
            f"degenerate n={n}: dI_TT win {res[f'n{n}']['dtt_win']:.2e} (law {res[f'n{n}']['dlaw_win']:.1e}), "
            f"out {res[f'n{n}']['dtt_out']:.2e} | other n={o}: dI_TT win {res[f'n{o}']['dtt_win']:.2e} (law {res[f'n{o}']['dlaw_win']:.2e}), "
            f"TT-law diff {res[f'n{o}']['d_tt_vs_law']:.1e} | law err n1/n2 {res['n1']['law_err']:.1e}/{res['n2']['law_err']:.1e} "
            f"| conv {res['conv'][0]:.1e}/{res['conv'][1]:.1e} thr {res['thresh']:.1e} -> {'BROKEN' if res['broken'] else 'holds'} ({res['time_s']:.0f} s)",
            flush=True,
        )
    return res


def section_run(S_sub=S_SCAN, explore=True):
    out = dict(S_scan=S_SCAN, thresh=THRESH, cases=[], exploratory=[])
    if OUT.exists():
        prev = json.load(open(OUT))
        out["cases"] = [r for r in prev.get("cases", []) if r["S"] not in S_sub]
    t0 = time.time()
    for S in S_sub:
        for key in ("a", "a2", "b"):
            if key == "a2" and S < 4:
                continue
            r = one_case(S, key)
            if r is None:
                print(f"S={S:3g} {key}: no admissible interval")
                continue
            out["cases"].append(r)
    out["cases"].sort(key=lambda r: (r["S"], r["key"]))
    broken = [(r["S"], r["key"]) for r in out["cases"] if r["broken"]]
    print("verdict:", f"BROKEN at {broken}" if broken else "class holds for all S <= 6, both orders")
    out["verdict"] = "broken" if broken else "holds"
    out["broken_at"] = broken
    if explore:
        print("exploratory (not under the rule): R = 0.1")
        out["exploratory"].append(one_case(3.0, "a", alpha=1.0, aac=10.0))
        out["exploratory"].append(one_case(6.0, "a2", alpha=1.0, aac=10.0))
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT.name} ({time.time() - t0:.0f} s)")


def section_runA():
    """First half (S <= 3), to keep each run under two minutes (BACKLOG rule 17a)."""
    section_run(S_sub=(1.5, 2.0, 3.0), explore=False)


def section_runB():
    """Second half (S >= 4) plus the exploratory rows."""
    section_run(S_sub=(4.0, 5.0, 6.0), explore=True)


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    ax = axs[0]
    for key, mk, lab in (("a", "o", "same $\\Theta_1$ (J$_1$ rising)"), ("a2", "s", "same $\\Theta_1$ (J$_1$ descending)"), ("b", "^", "same $\\Theta_2$ (J$_2$ rising)")):
        rows = [r for r in d["cases"] if r["key"] == key]
        if not rows:
            continue
        S = [r["S"] for r in rows]
        n = rows[0]["n"]
        ax.plot(S, [r[f"n{n}"]["dtt_win"] for r in rows], mk + "-", label=f"{lab}: TT $|\\Delta I_{n}|$ in window")
        ax.plot(S, [r[f"n{3 - n}"]["dtt_win"] for r in rows], mk + ":", mfc="none", color="0.5", label=f"{lab}: other order $|\\Delta I_{3 - n}|$")
    ax.axhline(d["thresh"], color="C3", lw=1, ls="--")
    ax.text(1.6, d["thresh"] * 1.2, "threshold 5e-4", color="C3", fontsize=8)
    ax.set(yscale="log", xlabel="S", ylabel="max difference between the two histories", title="degenerate order vs the other order, full TT", ylim=(1e-8, 1))
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(color="#e4e3df", lw=0.6)
    ax = axs[1]
    for key, mk in (("a", "o"), ("a2", "s"), ("b", "^")):
        rows = [r for r in d["cases"] if r["key"] == key]
        if not rows:
            continue
        S = [r["S"] for r in rows]
        ax.plot(S, [max(r["n1"]["law_err"], r["n2"]["law_err"]) for r in rows], mk + "-", label=f"|I_TT − I_law| ({key})")
        ax.plot(S, [max(r["conv"]) for r in rows], mk + ":", color="0.5", mfc="none", label=f"convergence ({key})")
    ax.set(yscale="log", xlabel="S", ylabel="max over t, both histories, both orders", title="law error and convergence error")
    ax.legend(fontsize=7)
    ax.grid(color="#e4e3df", lw=0.6)
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
