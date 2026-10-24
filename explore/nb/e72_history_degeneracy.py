"""Entry 72: which acoustic histories does time-resolved diffraction NOT distinguish? (directive B1)

Steady-state diffraction of a filled constant-amplitude wave knows only (S, q); the transient knows the
amplitude history A(s) through the accumulated coupling of Entries 28-33,
    Theta_n(t) = int_0^1 J_n(S A(t - zeta)) d zeta,      I_n = sin^2(pi alpha Theta_n)   (Laue),
i.e. order n reads only the running mean of J_n(S A) over the last transit. Two histories whose J_1(S A)
differ by a function with zero mean over every window of length t_tr are therefore predicted to be
indistinguishable at order 1 while the window contains the whole perturbation, and distinguishable at
order 2 (J_2 of the same amplitudes is not zero-mean). This is the B1 question made quantitative, and
a falsification test of claim #21 (Entry 33, "Theta_n exact") on a non-monotone envelope.

Construction: A_1(s) = clip(s, 0, 1) (linear ramp over t_tr, then held); A_2 defined by
    J_1(S A_2(s)) = J_1(S A_1(s)) + delta(s),   delta = eps sin(2 pi m (s - a)/(b - a)) on [a, b], else 0,
inverted on the monotone branch S A <= 1.84 (S = 1.5). With [a, b] = [0.2, 0.7] and integer m, the window
[t - 1, t] contains all of delta for t in [0.7, 1.2]: there Theta_1 is identical by construction.

Pre-registered rule (before the TT run). Full TT, Laue, alpha = 0.5, alpha_ac = 40 (R = 0.0125, the Entry
67 conditions), 8 carrier phases at fixed envelope time, eta = +q and +2q:
  (i) |I_1^(2) - I_1^(1)| < 1e-3 (10x the Entry 67 law error) for all t in [0.7, 1.2]  -> order-1 degeneracy
      confirmed; otherwise claim #21 fails on this envelope class;
  (ii) the TT difference traces agree with the law's difference traces to 1e-3 at every t, both orders;
  (iii) max_t |I_2^(2) - I_2^(1)| in [0.7, 1.2] is compared with the law value; if the law predicts > 3e-3
       and TT confirms, order 2 breaks the degeneracy at this S; if the law predicts < 3e-3 the question is
       "not resolved at this S" and no claim is made.
Usage: python e72_history_degeneracy.py [law] [tt] [fig]
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.core import required_nz, rocking  # noqa: E402

OUT = HERE / "e72_history_degeneracy.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e72_history_degeneracy.png"
ALPHA, AAC, S = 0.5, 40.0, 1.5
Q = 2 * np.pi * AAC
EPS, M, A_, B_ = 0.1, 1, 0.2, 0.7
TIMES = np.linspace(0.0, 1.5, 151)
WINDOW = (0.7, 1.2)
_x = np.linspace(0.0, 1.8411, 40001)
_j = jv(1, _x)


def inv_j1(y):
    return np.interp(y, _j, _x)


def A1(s):
    return np.clip(s, 0.0, 1.0)


def delta(s):
    s = np.asarray(s, float)
    inside = (s >= A_) & (s <= B_)
    return np.where(inside, EPS * np.sin(2 * np.pi * M * (s - A_) / (B_ - A_)), 0.0)


def A2(s):
    s = np.asarray(s, float)
    return np.where(s >= 0, inv_j1(jv(1, S * A1(s)) + delta(s)) / S, 0.0)


def theta(env, n, t, nz=8001):
    z = np.linspace(0.0, 1.0, nz)
    return np.array([np.trapezoid(jv(n, S * env(x - z)), z) for x in np.atleast_1d(t)])


def law(env, n, t):
    return np.sin(np.pi * ALPHA * theta(env, n, t)) ** 2


def tt_trace(env, n, times, phases=8, nz=None):
    nz = nz or required_nz(ALPHA, AAC, S=S, eta_max=2 * Q)
    tt = np.asarray(times, float)[:, None, None]
    pp = (np.arange(phases) * 2 * np.pi / phases)[None, :, None]

    def phi(z):
        return S * env(tt - z) * np.sin(Q * (z - tt) + pp)

    inten = rocking(ALPHA, np.array([n * Q]), phi=phi, n_z=nz)
    return inten[..., 0].mean(axis=1)


def section_law():
    win = (TIMES >= WINDOW[0]) & (TIMES <= WINDOW[1])
    out = {"t": TIMES.tolist()}
    for n in (1, 2):
        L1, L2 = law(A1, n, TIMES), law(A2, n, TIMES)
        d = L2 - L1
        out[f"law_n{n}"] = [L1.tolist(), L2.tolist()]
        print(
            f"law n={n}: max|dI| in window {np.max(np.abs(d[win])):.2e}, outside {np.max(np.abs(d[~win])):.2e};"
            f"  I range [{L1.min():.3f}, {L1.max():.3f}]"
        )
    s = np.linspace(0, 1.2, 7)
    print("envelopes at s =", s, "\n  A1", np.round(A1(s), 3), "\n  A2", np.round(A2(s), 3))
    print(
        f"max |A2 - A1| = {np.max(np.abs(A2(np.linspace(0, 1.5, 3001)) - A1(np.linspace(0, 1.5, 3001)))):.3f}"
    )
    json.dump(out, open(OUT, "w"))


def section_tt():
    out = json.load(open(OUT))
    win = (TIMES >= WINDOW[0]) & (TIMES <= WINDOW[1])
    nz = required_nz(ALPHA, AAC, S=S, eta_max=2 * Q)
    print(f"n_z = {nz}")
    for n in (1, 2):
        T1, T2 = tt_trace(A1, n, TIMES, nz=nz), tt_trace(A2, n, TIMES, nz=nz)
        L1, L2 = np.array(out[f"law_n{n}"])
        d_tt, d_law = T2 - T1, L2 - L1
        out[f"tt_n{n}"] = [T1.tolist(), T2.tolist()]
        print(
            f"TT n={n}: max|dI_TT| in window {np.max(np.abs(d_tt[win])):.2e} (law {np.max(np.abs(d_law[win])):.2e}),"
            f" outside {np.max(np.abs(d_tt[~win])):.2e} (law {np.max(np.abs(d_law[~win])):.2e});"
            f" max|dI_TT - dI_law| {np.max(np.abs(d_tt - d_law)):.2e};"
            f" max|I_TT - I_law| env1 {np.max(np.abs(T1 - L1)):.2e} env2 {np.max(np.abs(T2 - L2)):.2e}"
        )
    # convergence at three times: phases x2, n_z x2
    sub = np.array([0.5, 0.95, 1.4])
    for n in (1, 2):
        a = tt_trace(A2, n, sub, phases=8, nz=nz)
        b = tt_trace(A2, n, sub, phases=16, nz=2 * nz)
        print(f"convergence n={n} (phases 8->16, n_z x2): max change {np.max(np.abs(a - b)):.2e}")
        out[f"conv_n{n}"] = float(np.max(np.abs(a - b)))
    json.dump(out, open(OUT, "w"))


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    t = np.array(d["t"])
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
    s = np.linspace(0, 1.5, 1501)
    axs[0].plot(s, A1(s), "k-", label="history 1: ramp")
    axs[0].plot(s, A2(s), "C3-", label="history 2: ramp + zero-mean wobble in J_1")
    axs[0].set(xlabel="retarded time s / t_tr", ylabel="envelope A(s)", title="two drive histories")
    axs[0].legend(fontsize=8)
    for ax, n in zip(axs[1:], (1, 2), strict=True):
        L1, L2 = np.array(d[f"law_n{n}"])
        ax.plot(t, L1, "k-", lw=1, label="law, history 1")
        ax.plot(t, L2, "C3-", lw=1, label="law, history 2")
        if f"tt_n{n}" in d:
            T1, T2 = np.array(d[f"tt_n{n}"])
            ax.plot(t[::3], T1[::3], "ko", ms=3, mfc="none", label="TT, history 1")
            ax.plot(t[::3], T2[::3], "C3s", ms=3, mfc="none", label="TT, history 2")
        ax.axvspan(*WINDOW, color="#e4e3df", lw=0, label="window holds the whole wobble")
        ax.set(xlabel="t / t_tr", ylabel=f"I_{n}(t), phase-averaged", title=f"order n = {n}")
        ax.legend(fontsize=7)
    for ax in axs:
        ax.grid(color="#e4e3df", lw=0.6)
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["law", "tt", "fig"]:
        globals()[f"section_{sec}"]()
