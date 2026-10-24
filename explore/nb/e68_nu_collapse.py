"""Entry 68: is the linearity limit of the Fourier-weight estimator set by one number?

Central-claim test (notes/paper_candidates.md). The claim: the satellite's own dynamical phase
    nu = first-order satellite amplitude at the surface = max_t |R_lin(t) - r0|
decides whether time-resolved Bragg reflectivity is a linear spectrometer of the strain (nu << 1,
Entry 59/65 transfer function) or a dynamical satellite (nu ~ 1, Entry 1/9 formation law), and it does
so regardless of how the strained depth is produced (packet length vs penetration depth) and of the
absorption. Entry 65 had this on four points. Here it is a scan.

Setup: semi-infinite symmetric Bragg crystal, the retarded-time ODE of Entry 66 (dR/ds, one integration
gives the whole R(t)); length unit Lambda (alpha = 1), time unit Lambda / v, absorption kappa = 2 pi g.
Packet phi = S g(s), g = exp(-s / tau_d) sin(q0 s), tau_d = ncyc / x, q0 = 2 pi x, x = satellite position
in Darwin half-widths. Measure at eta = q0 (the spectral peak) and at 0.8 q0, 1.25 q0:
    K(S) = (W(S) / S) / (W(S_lin) / S_lin),   W = |Fourier component of |R|^2 at omega = eta|.
K = 1 is linear. The kinematic Bessel factor 2 J_1(S) / S is reported separately so that the dynamical
part (K / kinematic) can be read off; for S <= 0.3 it is within 1.1 % of 1 anyway.

Falsification rule, fixed before the run: if at equal nu the quantity 1 - K spreads by more than a factor 2
across (g, x, ncyc), the single-parameter claim is false in its strong form.

Usage: python e68_nu_collapse.py [run] [fig] [conv]
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e65_si400 import root  # noqa: E402

OUT = HERE / "e68_nu_collapse.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e68_nu_collapse.png"
S_LIN = 1e-3
S_GRID = [S_LIN, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0]
G_GRID = [0.01, 0.046, 0.12]
X_GRID = [2.0, 5.0]
NCYC_GRID = [1.0, 3.0, 10.0, 30.0]
ETA_FRAC = [0.8, 1.0, 1.25]


def surface_R_multi(kappa, eta, S, g, T, per_unit):
    """R(t) at the surface, vectorised over S (rows) and eta (columns); alpha = 1."""
    eta = np.asarray(eta, float)
    S = np.asarray(S, float)[:, None]
    c = 1j * np.pi
    ie = 1j * (eta + 1j * kappa)[None, :]
    R = np.broadcast_to(
        np.array([root(1.0, e, kappa) for e in eta])[None, :], (S.size, eta.size)
    ).copy()
    n = int(np.ceil(T * per_unit))
    h = T / n
    s = np.arange(n + 1) * h

    def f(si, r):
        e = np.exp(1j * S * g(si))
        return c * np.conj(e) + (ie + c * e * r) * r

    out = np.empty((n + 1, S.size, eta.size), complex)
    out[0] = R
    for i in range(n):
        k1 = f(s[i], R)
        k2 = f(s[i] + h / 2, R + h / 2 * k1)
        k3 = f(s[i] + h / 2, R + h / 2 * k2)
        k4 = f(s[i] + h, R + h * k3)
        R = R + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = R
    return s, out


def weight(t, R, eta):
    inten = np.abs(R) ** 2
    inten = inten - inten.mean(axis=0)
    w = np.full(t.size, t[1] - t[0])
    w[[0, -1]] /= 2
    ph = np.exp(1j * np.abs(eta)[None, None, :] * t[:, None, None])
    return np.abs(np.sum(w[:, None, None] * inten * ph, axis=0))


def packet(q0, ncyc):
    tau_d = ncyc * 2 * np.pi / q0
    return lambda s: np.where(
        s > 0, np.exp(-np.clip(s, 0, None) / tau_d) * np.sin(q0 * s), 0.0
    ), tau_d


def penetration(x, kappa):
    """1/e amplitude depth of the first-order satellite source: Re of Lambda - i W (Entry 59)."""
    W = 2 * np.pi * x
    R0 = root(1.0, W, kappa)
    lam = 1j * (W + 1j * kappa) + 2j * np.pi * R0
    return 1.0 / abs((lam - 1j * W).real)


def case(g, x, ncyc, per_cycle=64, S_grid=S_GRID, T_scale=1.0):
    kappa = 2 * np.pi * g
    q0 = 2 * np.pi * x
    gfun, tau_d = packet(q0, ncyc)
    ell_p = penetration(x, kappa)
    eta = q0 * np.array(ETA_FRAC)
    T = (5 * tau_d + 4 * ell_p + 1.0) * T_scale
    per_unit = per_cycle * (eta.max() + max(S_grid) * q0) / (2 * np.pi)
    t, R = surface_R_multi(kappa, eta, S_grid, gfun, T, per_unit)
    W = weight(t, R, eta)  # (nS, n_eta)
    r0 = R[0]
    nu_lin = np.max(np.abs(R[:, 0, :] - r0[0]), axis=0)  # per eta, at S_LIN
    K = (W / np.asarray(S_grid)[:, None]) / (W[0] / S_LIN)[None, :]
    nu = np.asarray(S_grid)[:, None] / S_LIN * nu_lin[None, :]
    return dict(
        g=g,
        x=x,
        ncyc=ncyc,
        kappa=kappa,
        q0=q0,
        tau_d=tau_d,
        ell_p=ell_p,
        T=T,
        n_steps=int(t.size - 1),
        eta_frac=ETA_FRAC,
        S=list(S_grid),
        K=K.tolist(),
        nu=nu.tolist(),
        W_lin_over_S=(W[0] / S_LIN).tolist(),
        kin=(2 * jv(1, np.asarray(S_grid)) / np.asarray(S_grid)).tolist(),
    )


def section_run():
    rows = []
    for g in G_GRID:
        for x in X_GRID:
            for ncyc in NCYC_GRID:
                r = case(g, x, ncyc)
                rows.append(r)
                K = np.array(r["K"])[:, 1]
                nu = np.array(r["nu"])[:, 1]
                print(
                    f"g={g:<6} x={x:<3} ncyc={ncyc:<5} ell_p={r['ell_p']:.2f} tau_d={r['tau_d']:.2f} "
                    f"steps={r['n_steps']}: "
                    + " ".join(
                        f"S={s:g}:nu={n:.3f},K={k:.3f}"
                        for s, n, k in zip(r["S"][1:], nu[1:], K[1:], strict=True)
                    ),
                    flush=True,
                )
    json.dump(dict(S_lin=S_LIN, rows=rows), open(OUT, "w"), indent=1)


def section_conv():
    """Step and window convergence at one strongly nonlinear case."""
    for per in (64, 128):
        r = case(0.01, 2.0, 30.0, per_cycle=per, S_grid=[S_LIN, 0.1, 0.3])
        print(f"per_cycle={per}: K = {np.array(r['K'])[1:, 1]}  nu = {np.array(r['nu'])[1:, 1]}")


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    rows = d["rows"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    colors = {0.01: "#2a78d6", 0.046: "#eb6834", 0.12: "#1baf7a"}
    markers = {1.0: "o", 3.0: "s", 10.0: "^", 30.0: "D", 100.0: "x"}
    nu_ref = np.geomspace(1e-2, 3, 200)
    for ax, sel in zip(axs, (slice(1, 7), slice(1, None)), strict=True):
        for r in rows:
            K = np.array(r["K"])[sel, 1]
            nu = np.array(r["nu"])[sel, 1]
            ax.plot(
                nu,
                1 - K,
                markers[r["ncyc"]],
                color=colors[r["g"]],
                ms=5,
                mfc="none" if r["x"] == 5 else None,
                lw=0,
            )
        ax.plot(
            nu_ref,
            1 - np.tanh(nu_ref) ** 2 / nu_ref**2,
            "k-",
            lw=1,
            label=r"Bragg law $1-\tanh^2\nu/\nu^2$",
        )
        ax.plot(
            nu_ref,
            1 - np.sin(nu_ref) ** 2 / nu_ref**2,
            "k--",
            lw=1,
            label=r"Laue law $1-\sin^2\nu/\nu^2$",
        )
        ax.set(
            xscale="log",
            yscale="log",
            xlabel=r"$\nu$ = first-order satellite amplitude at the surface",
            ylabel="1 - K  (saturation of the weight at the spectral peak)",
            ylim=(1e-3, 2),
        )
        ax.grid(color="#e4e3df", lw=0.6)
    axs[0].set_title("S <= 0.3 (kinematic Bessel factor < 1.1 %)", fontsize=10)
    axs[1].set_title("all S (includes S = 0.5, 1: Bessel saturation adds)", fontsize=10)
    from matplotlib.lines import Line2D

    h = [Line2D([], [], color=c, marker="o", lw=0, label=f"g = {g}") for g, c in colors.items()]
    h += [
        Line2D([], [], color="k", marker=m, lw=0, label=f"{n:g} cycles") for n, m in markers.items()
    ]
    h += [
        Line2D([], [], color="k", marker="o", lw=0, mfc="none", label="open: x = 5; filled: x = 2")
    ]
    axs[1].legend(
        handles=h + axs[1].get_legend_handles_labels()[0], fontsize=7, loc="lower right", ncol=2
    )
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


def section_stats():
    """Spread of the saturation at fixed nu across (g, x, ncyc): the falsification rule of the docstring."""
    d = json.load(open(OUT))
    for lab, corr in (("1 - K", False), ("1 - K / (2 J_1(S) / S)", True)):
        print(f"\n{lab} interpolated at fixed nu (log-log within each series, S <= 1):")
        print("  nu    min    max   max/min   n_series  (least / most saturating series)")
        for nu0 in (0.1, 0.2, 0.3, 0.5, 0.8):
            vals = []
            for r in d["rows"]:
                nu = np.array(r["nu"])[1:, 1]
                K = np.array(r["K"])[1:, 1]
                if corr:
                    K = K / np.array(r["kin"])[1:]
                y = 1 - K
                ok = (y > 0) & (nu > 0)
                if ok.sum() >= 2 and nu[ok].min() <= nu0 <= nu[ok].max():
                    v = np.exp(np.interp(np.log(nu0), np.log(nu[ok]), np.log(y[ok])))
                    vals.append((v, f"g={r['g']},x={r['x']:g},n={r['ncyc']:g}"))
            if len(vals) < 2:
                continue
            vals.sort()
            print(
                f"  {nu0:<5} {vals[0][0]:.4f} {vals[-1][0]:.4f} {vals[-1][0] / vals[0][0]:6.2f}   {len(vals):2d}"
                f"      {vals[0][1]} / {vals[-1][1]}"
            )
    print("\nSame nu, different route (S <= 0.3 only): short vs long packet at g = 0.01, x = 5")
    for r in d["rows"]:
        if r["g"] == 0.01 and r["x"] == 5.0:
            nu = np.array(r["nu"])[1:7, 1]
            K = np.array(r["K"])[1:7, 1]
            print(
                f"  ncyc={r['ncyc']:>4g}: "
                + " ".join(f"({n:.2f},{1 - k:.3f})" for n, k in zip(nu, K, strict=True))
            )


def section_fit2():
    """Constructive follow-up to the failed collapse: is 1 - K = a S^2 + b nu^2 (local + accumulated)?

    Least squares over S <= 0.3 at the spectral peak, then the same with b fixed to the Bragg law 2/3.
    """
    d = json.load(open(OUT))
    S, nu, y, lab = [], [], [], []
    for r in d["rows"]:
        for i in range(1, 7):
            S.append(r["S"][i])
            nu.append(r["nu"][i][1])
            y.append(1 - r["K"][i][1])
            lab.append(f"g={r['g']},x={r['x']:g},n={r['ncyc']:g},S={r['S'][i]}")
    S, nu, y = map(np.array, (S, nu, y))
    A = np.stack([S**2, nu**2], axis=1)
    (a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ [a, b]
    rel = pred / y
    print(f"fit 1 - K = a S^2 + b nu^2 over {y.size} points (S <= 0.3): a = {a:.4f}, b = {b:.4f}")
    print(
        f"  pred / actual: median {np.median(rel):.3f}, 5-95 % [{np.percentile(rel, 5):.3f}, {np.percentile(rel, 95):.3f}], "
        f"min {rel.min():.3f} ({lab[rel.argmin()]}), max {rel.max():.3f} ({lab[rel.argmax()]})"
    )
    big = y > 0.01
    print(
        f"  where 1 - K > 0.01 ({big.sum()} points): pred / actual min {rel[big].min():.3f}, max {rel[big].max():.3f}"
    )
    for bb in (2 / 3,):
        aa = np.sum(S**2 * (y - bb * nu**2)) / np.sum(S**4)
        rel2 = (aa * S**2 + bb * nu**2) / y
        print(
            f"  with b = {bb:.3f} fixed (Bragg tanh law): a = {aa:.4f}; pred / actual 5-95 % "
            f"[{np.percentile(rel2, 5):.3f}, {np.percentile(rel2, 95):.3f}], where 1 - K > 0.01: "
            f"[{rel2[big].min():.3f}, {rel2[big].max():.3f}]"
        )
    nu_only = np.sum(nu**2 * y) / np.sum(nu**4)
    rel3 = nu_only * nu**2 / y
    print(
        f"  nu-only fit b = {nu_only:.3f}: pred / actual 5-95 % [{np.percentile(rel3, 5):.3f}, {np.percentile(rel3, 95):.3f}]"
    )


def section_window():
    """Window convergence: the weight window doubled at a short and a long packet (g = 0.01, x = 5)."""
    for ncyc in (1.0, 30.0):
        for ts in (1.0, 2.0):
            r = case(0.01, 5.0, ncyc, S_grid=[S_LIN, 0.1, 0.3], T_scale=ts)
            print(
                f"ncyc={ncyc:g} T={r['T']:.1f}: K = {np.array(r['K'])[1:, 1]}  nu = {np.array(r['nu'])[1:, 1]}"
            )


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
