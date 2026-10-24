"""Entry 83 (T2m): the beam condition of Entry 82 at GHz-train frequencies, x = f / f_D = 1-10.

Same machinery as e82_divergence.py (full TT beam average of the Entry 66 packet at S/1000 against the convolved
closed form CF-B), at Si 400 / 17.48 keV and 8 keV, x = 1, 2, 5, 10, windows 52 and 300 ns, Gaussian beam spread
sigma = 0.03-2 Darwin half-widths. Pre-registered rule (notebook Entry 83): for each (E, x, T_w) the largest sigma at
which |H_conv / H_single - 1| <= 5 % (interpolated), in Darwin half-widths, urad and meV; "GHz train measurable with an
ordinary beam" if that exceeds a few urad / a few meV at x = 1-5, else "single-ray condition stands".
Usage: python e83_lowx.py [run] [fig]
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from e65_si400 import SI400, params  # noqa: E402
from e75_window_cf import integrate_packet, setup  # noqa: E402
from e82_divergence import cf_complex, delta_grid, fourier_complex, norm  # noqa: E402

OUT = HERE / "e83_lowx.json"
FIG = HERE.parents[1] / "figures" / "paper" / "e83_lowx.png"
XS = np.array([1.0, 2.0, 5.0, 10.0])
SIGMAS = (0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0)
T_NS = (52.0, 300.0)
TOL = 0.05


def angle_scale(E):
    th, xr, _ = SI400[E]
    thD = xr / np.sin(np.radians(2 * th))  # Darwin half-width, rad
    return thD, thD / np.tan(np.radians(th)) * E * 1e6  # urad per half-width -> (rad, meV per half-width)


def section_run():
    t0 = time.time()
    out = dict(xs=XS.tolist(), sigmas=SIGMAS, T_ns=T_NS, tol=TOL, cases={}, worst_tt_cf=0.0)
    for E in (17.48, 8.0):
        base = setup(E, s_lin=1000.0)
        kappa = base["kappa"]
        g = params(E)["g"]
        thD, meV = angle_scale(E)
        eta_nom = 2 * np.pi * XS
        f_GHz = XS * params(E)["f_D_GHz"]
        # single ray
        st0 = dict(base, eta=eta_nom)
        run0 = integrate_packet(st0)
        H1 = {}
        for T_ns in T_NS:
            T = T_ns * 1e-9 / base["t_unit"]
            I0, Ta = fourier_complex(st0, run0, T, eta_nom)
            H1[T_ns] = np.abs(I0) / norm(base, eta_nom, Ta)
        rows = {}
        for sigma in SIGMAS:
            sig_eta = 2 * np.pi * sigma
            d, p = delta_grid(sig_eta, kappa)
            eta_all = (eta_nom[:, None] + d[None, :]).ravel()
            om_all = np.repeat(eta_nom, d.size)
            st = dict(base, eta=eta_all)
            run = integrate_packet(st)
            rows[sigma] = {}
            for T_ns in T_NS:
                T = T_ns * 1e-9 / base["t_unit"]
                Itt, Ta = fourier_complex(st, run, T, om_all)
                Icf = cf_complex(st, eta_all, om_all, Ta)
                Wtt = np.abs((Itt.reshape(XS.size, -1) * p[None, :]).sum(axis=1))
                Wcf = np.abs((Icf.reshape(XS.size, -1) * p[None, :]).sum(axis=1))
                out["worst_tt_cf"] = max(out["worst_tt_cf"], float(np.max(np.abs(Wtt / Wcf - 1))))
                rows[sigma][T_ns] = (Wtt / norm(base, eta_nom, Ta) / H1[T_ns]).tolist()
        # largest sigma within 5 %, per (x, T_w), by linear interpolation on the scan
        smax = {}
        for T_ns in T_NS:
            smax[T_ns] = []
            for i in range(XS.size):
                r = np.array([abs(rows[s][T_ns][i] - 1) for s in SIGMAS])
                sg = np.array(SIGMAS)
                if r[0] > TOL:
                    smax[T_ns].append(0.0)
                elif r[-1] <= TOL:
                    smax[T_ns].append(float(sg[-1]))
                else:
                    j = int(np.argmax(r > TOL))
                    smax[T_ns].append(float(np.interp(TOL, [r[j - 1], r[j]], [sg[j - 1], sg[j]])))
        out["cases"][f"E{E:g}"] = dict(g=g, f_GHz=f_GHz.tolist(), urad_per_hw=thD * 1e6, meV_per_hw=meV,
                                        H_single={str(k): v.tolist() for k, v in H1.items()},
                                        ratio={f"s{s:g}": {str(k): v for k, v in rows[s].items()} for s in SIGMAS},
                                        sigma_max_hw={str(k): v for k, v in smax.items()},
                                        expect_hw=(0.25 * (g + 1 / (2 * XS))).tolist())
        print(f"Si 400 {E:g} keV: g = {g:.4f}, Darwin half-width {thD * 1e6:.2f} urad = {meV:.1f} meV; f = " + " ".join(f"{v:.2f}" for v in f_GHz) + " GHz for x = 1, 2, 5, 10")
        for T_ns in T_NS:
            print(f"  T_w = {T_ns:g} ns: H_single = " + " ".join(f"{v:.3f}" for v in H1[T_ns]))
            for s in SIGMAS:
                print(f"     sigma = {s:4g} hw: H_conv/H_single = " + " ".join(f"{v:.3f}" for v in rows[s][T_ns]))
            sm = np.array(smax[T_ns])
            print(f"     sigma_max (5 %) = " + " ".join(f"x={x:g}:{v:.3f} hw ({v * thD * 1e6:.2f} urad, {v * meV:.1f} meV)" for x, v in zip(XS, sm, strict=True))
                  + "   expected 0.25 (g + 1/2x) = " + " ".join(f"{v:.3f}" for v in 0.25 * (g + 1 / (2 * XS))))
        print(f"  ({time.time() - t0:.0f} s)", flush=True)
    # convergence at 8 keV, x = 2, sigma = 0.3, 300 ns
    base = setup(8.0, s_lin=1000.0)
    for du in (0.025, 0.0125):
        d, p = delta_grid(2 * np.pi * 0.3, base["kappa"], du=du)
        st = dict(base, eta=2 * np.pi * 2.0 + d)
        run = integrate_packet(st)
        Itt, _ = fourier_complex(st, run, 300e-9 / base["t_unit"], np.full(d.size, 2 * np.pi * 2.0))
        W = abs((Itt * p).sum())
        print(f"convergence 8 keV, x = 2, sigma = 0.3, 300 ns, du = {du}: {d.size} rays, W_conv = {W:.6e}")
        out[f"conv_du{du}"] = float(W)
    print(f"TT vs closed-form convolution, worst: {out['worst_tt_cf'] * 100:.3f} %")
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"wrote {OUT.name} ({time.time() - t0:.0f} s)")


def section_fig():
    import matplotlib.pyplot as plt

    d = json.load(open(OUT))
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, E in zip(axs, ("17.48", "8"), strict=True):
        c = d["cases"][f"E{E}"]
        for i, (x, col) in enumerate(zip(d["xs"], ("C0", "C1", "C3", "k"), strict=True)):
            ax.plot(d["sigmas"], [c["ratio"][f"s{s:g}"]["300.0"][i] for s in d["sigmas"]], "o-", color=col, label=f"x = {x:g} ({c['f_GHz'][i]:.1f} GHz)")
        ax.axhspan(1 - d["tol"], 1 + d["tol"], color="#e4e3df", lw=0, zorder=0)
        ax.set(xscale="log", xlabel="beam spread $\\sigma$ (Darwin half-widths)", ylabel="$H_{conv}/H_{single}$ (300 ns)", title=f"Si 400, {E} keV: {c['urad_per_hw']:.1f} µrad, {c['meV_per_hw']:.0f} meV per half-width", ylim=(0, 1.1))
        ax.legend(fontsize=8, frameon=False)
        ax.grid(color="#e4e3df", lw=0.6)
    fig.savefig(FIG, dpi=170)
    print("wrote", FIG)


if __name__ == "__main__":
    for sec in sys.argv[1:] or ["run", "fig"]:
        globals()[f"section_{sec}"]()
