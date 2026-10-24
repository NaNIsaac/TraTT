"""Paper figures that are not yet files (BACKLOG T2i; notes/paper_draft.md section 3). No new computation:
every panel is drawn from an existing json / npz, or from the 0.2 s packet run of e75 (Fig. 1a, R(t) itself).
Poster figures stay frozen; outputs go to figures/paper/.

    fig1  (a) R(t) of the Entry 66 packet at one detuning (e75 run, Si 400 / 17.48 keV, S/100) with the
              Fourier component at omega = |eta| marked; (b) W / kinematic law against omega / omega_D for
              three packet spectra at kappa / 2 pi alpha = 0.12 (e58_abs.json, e58_darwin.json; Entry 59).
    fig2  |F_pred| / |F_TT| against eta / (2 pi alpha) for the plane-wave S&W form, its normalised form and the
          Bloch-wave (Riccati) form, R_D = 0.05 and 0.2 (e64_sw.json; Entries 64, 69).
    fig4b the formation-law collapse I_n(alpha |J_n(S)| t) for Laue (sin^2) and Bragg (tanh^2), four cases
          (poster_f1_points.npz; Entries 1, 9), re-drawn in the paper style.
Captions with data file, script and commit: figures/paper/CAPTIONS.md (written by this script).
Usage: python paper_figs.py [fig1] [fig2] [fig4b]
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTD = ROOT / "figures" / "paper"
sys.path.insert(0, str(HERE))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GRID = dict(color="#e4e3df", lw=0.6)
COMMIT = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
CAPTIONS = {}


def fig1():
    from e59_transfer import G, packet, predict
    from e75_window_cf import integrate_packet, setup

    # (a) R(t) at 40 GHz detuning, 17.48 keV, S/100; window 1 ns; the Fourier component at omega = |eta|
    st = setup(17.48)
    f_sel = 40.0
    j = list(np.round(st["eta"] / (2 * np.pi * 1e9 * st["t_unit"]))).index(f_sel)
    run = integrate_packet(st)
    s, R, r0 = run["s"], run["R"][:, j], run["r1"][j]
    t_ps = s * st["t_unit"] * 1e12
    dI = np.abs(R) ** 2 - abs(r0) ** 2
    eta = st["eta"][j]
    # continuation past the packet for the picture only (first order): delta ~ e^{Lambda (s - s_p)}
    s2 = np.linspace(s[-1], 1e-9 / st["t_unit"], 20000)
    u = (R[-1] - r0) / (R[-1] - run["r2"][j]) * np.exp(run["Lam"][j] * (s2 - s[-1]))
    R2 = (r0 - u * run["r2"][j]) / (1 - u)
    dI2 = np.abs(R2) ** 2 - abs(r0) ** 2
    tt = np.concatenate([t_ps, s2 * st["t_unit"] * 1e12])
    yy = np.concatenate([dI, dI2])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ax[0].plot(tt, yy / abs(r0) ** 2, color="#2a78d6", lw=1.0)
    ax[0].plot(tt, np.abs(np.cos(eta * np.concatenate([s, s2]))) * 0 + 0, lw=0)  # keeps limits
    ax[0].set(xlabel="time after the packet enters (ps)", ylabel=r"$(|R|^2 - |R_0|^2)\,/\,|R_0|^2$",
              title=f"a  Si 400, 17.48 keV, detuning {f_sel:g} GHz, peak strain 1e-5", xlim=(0, 1000))
    ax[0].text(0.98, 0.95, r"oscillation at $\omega = |\eta|$ = 40 GHz" + "\n" + r"amplitude decays as $e^{\mathrm{Re}\,\Lambda\, t}$ (tens of ns)",
               transform=ax[0].transAxes, ha="right", va="top", fontsize=9, color="#52514e")
    ax[0].grid(**GRID)
    # (b) the Entry 59 one-curve collapse, as in e59_transfer.figure panel b
    alpha, kappa, S = 8.0, 6.0, 0.02
    T = 0.6
    wd = 2 * np.pi * alpha
    D = json.load(open(HERE / "e58_darwin.json"))
    A = json.load(open(HERE / "e58_abs.json"))
    cases = [(0.05, A["S0.02"], 2 * np.pi * alpha / 0.05), (0.2, D["RD0.2"], D["RD0.2"]["q0"]), (0.5, D["RD0.5"], D["RD0.5"]["q0"])]
    for (RD, r, q0), c in zip(cases, ("#2a78d6", "#eb6834", "#1baf7a"), strict=True):
        e = np.array(r["eta"])
        Wn = np.array(r["W"])
        Wp, R0 = predict(alpha, kappa, S, q0, e)
        Wk = np.abs(R0) * np.abs(G(*packet(S, q0), e)) * np.pi * alpha * S / kappa / T
        x = e / wd
        ax[1].plot(x, Wn / Wk, "o", ms=3.5, color=c, mec="white", mew=0.4, label=rf"TT, $R_D$ = {RD}")
        ax[1].plot(x, Wp / Wk, "-", color=c, lw=1.4)
    ax[1].axhline(1, color="#52514e", lw=0.8, ls=":")
    ax[1].set(xscale="log", yscale="log", xlabel=r"$\omega/\omega_D$ = phonon frequency / Darwin half-width",
              ylabel="Fourier weight / kinematic law", title=r"b  one curve at $\kappa/2\pi\alpha$ = 0.12 (points TT, lines closed form)")
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].grid(**GRID)
    out = OUTD / "fig1_estimator.png"
    fig.savefig(out, dpi=200)
    CAPTIONS["fig1_estimator.png"] = (
        "Fig. 1. (a) Reflectivity change of a semi-infinite Si 400 crystal at 17.48 keV after the Entry 66 strain packet "
        "(87 GHz, 0.9 cycles, peak strain 1e-5) enters at t = 0, at a detuning of 40 GHz: the signal oscillates at the "
        "detuning frequency and its amplitude decays on the absorption/extinction scale (Entry 75; "
        "explore/nb/e75_window_cf.py integrate_packet, continuation past the packet at first order for the picture). "
        "(b) Fourier weight of R(t) at omega = |eta| divided by the kinematic law, against omega / omega_D, for three packet "
        "spectra (R_D = 0.05, 0.2, 0.5) at kappa / 2 pi alpha = 0.12: full TT (points, explore/nb/e58_abs.json, "
        "e58_darwin.json) and the first-order closed form (lines, explore/nb/e59_transfer.py predict); Entry 59. "
        f"Drawn by explore/nb/paper_figs.py fig1 at {COMMIT}."
    )
    print("wrote", out)


def fig2():
    d = json.load(open(HERE / "e64_sw.json"))
    WD = 2 * np.pi * 8.0
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True, sharey=True)
    for a, RD in zip(ax, (0.2, 0.05), strict=True):
        r = d[f"RD{RD}"]
        eta = np.array(r["eta"])
        Ft = np.array(r["TT"][0]) + 1j * np.array(r["TT"][1])
        big = np.abs(Ft) > 0.02 * np.abs(Ft).max()
        for kind, lab, c in (("SW", "plane-wave first order (S&W eq. 29)", "#eb6834"), ("SWn", "same, normalised to the incident wave", "#e3b505"), ("Riccati", "Bloch-wave (Riccati) first order, Entry 59", "#2a78d6")):
            Fp = np.array(r[kind][0]) + 1j * np.array(r[kind][1])
            a.plot(eta[big] / WD, np.abs(Fp[big]) / np.abs(Ft[big]), ".", ms=4, color=c, label=lab)
        a.axhline(1, color="#52514e", lw=0.8, ls=":")
        for x in (r["W"], -r["W"]):
            a.axvline(x / WD, color="0.7", lw=0.8, ls="--")
        a.set(xlabel=r"$\eta/(2\pi\alpha)$ (detuning in Darwin half-widths)", yscale="log",
              title=rf"$R_D$ = {RD}: satellite at $\pm${1 / RD:.0f} half-widths (dashed)")
        a.grid(**GRID)
    ax[0].set_ylabel(r"$|F_{\rm pred}|\,/\,|F_{\rm TT}|$ (Fourier amplitude of $R$ at $\omega = W$)")
    ax[0].legend(frameon=False, fontsize=8, loc="upper left")
    out = OUTD / "fig2_planewave_vs_bloch.png"
    fig.savefig(out, dpi=200)
    CAPTIONS["fig2_planewave_vs_bloch.png"] = (
        "Fig. 2. First-order predictions of the Fourier amplitude of R(t) at the mode frequency W, divided by the full TT value, "
        "against the detuning in Darwin half-widths (alpha = 8, kappa = 6, S = 0.02, single mode): the plane-wave form of "
        "Sondhauss and Wark (eq. 29), its incident-wave normalised form, and the Bloch-wave (Riccati) form of Entry 59. "
        "Left: mode at 5 half-widths (R_D = 0.2); right: 20 half-widths (R_D = 0.05). Only points with |F_TT| above 2 % of "
        "its maximum are shown. Data explore/nb/e64_sw.json from explore/nb/e64_sw_compare.py; Entries 64, 69 "
        f"(S&W transcription, PDF not seen). Drawn by explore/nb/paper_figs.py fig2 at {COMMIT}."
    )
    print("wrote", out)


def fig4b():
    from scipy.special import jv

    z = np.load(HERE / "poster_f1_points.npz")
    t, cases, laue, bragg = z["t"], z["cases"], z["laue"], z["bragg"]
    s = np.linspace(0, 2.3, 500)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True, sharey=True)
    for a, data, law, title in ((ax[0], laue, np.sin(np.pi * s) ** 2, r"Laue: $I_n = \sin^2(\pi\alpha|J_n(S)|\tilde t)$"),
                                (ax[1], bragg, np.tanh(np.pi * s) ** 2, r"Bragg: $I_n = \tanh^2(\pi\alpha|J_n(S)|\tilde t)$")):
        a.plot(s, law, color="#1d1d1b", lw=1.8, zorder=3, label="law, no fitted parameter")
        for (S, al, n), I, c, m in zip(cases, data, ("#2a78d6", "#eb6834", "#1baf7a", "#e3b505"), ("o", "s", "^", "D"), strict=True):
            a.plot(al * abs(jv(int(n), S)) * t, I, m, color=c, mec="white", mew=0.6, ms=6, ls="none", zorder=2,
                   label=rf"TT: $S$ = {S:g}, $\alpha$ = {al:g}, $n$ = {int(n)}")
        a.set(xlabel=r"$\alpha\,|J_n(S)|\,\tilde t$", title=title, xlim=(0, 2.3))
        a.grid(**GRID)
    ax[0].set_ylabel(r"satellite centre intensity $I_n$")
    ax[0].legend(frameon=False, fontsize=8)
    out = OUTD / "fig4b_formation_collapse.png"
    fig.savefig(out, dpi=200)
    CAPTIONS["fig4b_formation_collapse.png"] = (
        "Fig. 4b. Formation of the n-th acoustic satellite behind a step front: centre intensity against "
        "alpha |J_n(S)| t/t_tr for four (S, alpha, n) cases in Laue (left) and Bragg (right) geometry, full TT (markers) "
        "and the independent-satellite laws sin^2 and tanh^2 (line), no fitted parameter. Data explore/nb/poster_f1_points.npz "
        "(the cached TT points of the poster figure f1, explore/make_poster_figs.py; Entries 1, 9); the poster figure itself "
        f"is unchanged. Drawn by explore/nb/paper_figs.py fig4b at {COMMIT}."
    )
    print("wrote", out)


def write_captions():
    p = OUTD / "CAPTIONS.md"
    old = {}
    if p.exists():
        cur = None
        for line in p.read_text(encoding="utf-8").splitlines():
            if line.startswith("## "):
                cur = line[3:].strip()
                old[cur] = ""
            elif cur:
                old[cur] += line + "\n"
    old.update({k: v + "\n" for k, v in CAPTIONS.items()})
    text = "# Captions for figures/paper (data file, script, commit per figure)\n\n" + "".join(f"## {k}\n{v.strip()}\n\n" for k, v in old.items())
    p.write_text(text, encoding="utf-8", newline="\n")
    print("wrote", p)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["fig1", "fig2", "fig4b"]:
        globals()[sec]()
    write_captions()
