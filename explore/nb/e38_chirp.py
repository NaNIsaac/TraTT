"""Entry 38: chirped acoustic wave.

A1  spatial chirp q(zeta) = q0 [1 + beta (zeta - 1/2)], fully filled crystal (t = 1.5).
    Satellite n = 1 local resonance eta_c(zeta) = q(zeta) sweeps q0 beta across the slab.
    Landau-Zener prediction (two-beam, linear detuning): diabatic D0 -> Dh transfer
    P = 1 - exp(-2 pi a^2 / kappa), a = pi alpha |J_1(S)|, kappa = q0 beta, giving a flat
    plateau of height P over the swept eta range.  Kinematic limit: P ~ 2 pi a^2 / kappa.
A2  retarded-time chirp (frequency chirp at the transducer): transient trajectory of the
    satellite while the front fills the crystal.
E   "focusing": search a quadratic and a cubic phase for peak/integral and FWHM gains over
    the uniform-q control, TT and kinematic.
Usage: python e38_chirp.py [A1] [A2] [E]
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from qmod import PhaseDrive, features, kinematic, retarded_chirp, spatial_chirp, uniform  # noqa: E402
from tratt.core import required_nz, rocking  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

AAC, S = 40.0, 1.5
Q0 = 2 * np.pi * AAC
J1 = abs(jv(1, S))


def line(alpha, d, t, half=0.45, step=0.1, kin=True, nzmul=1, qmax=None):
    eta = np.arange(Q0 * (1 - half), Q0 * (1 + half), step)
    nz = nzmul * required_nz(alpha, (qmax or Q0) / (2 * np.pi), S=S, eta_max=eta.max())
    I = rocking(alpha, eta, phi=d.phi(t), n_z=nz)
    I0 = rocking(alpha, eta, n_z=nz)
    K = kinematic(alpha, d, eta, t) if kin else None
    return eta - Q0, I - I0, (K - kinematic(alpha, PhaseDrive(S=0.0, alpha_ac=AAC), eta, t)) if kin else None


def section_A1():
    print("== A1 spatial chirp, t = 1.5 (filled), n = 1, S = 1.5 ==")
    print(" alpha  beta | Gamma=a^2/k | LZ P   | TT plateau  peak  width  integral  npk | "
          "kin plateau  integral | uniform-q integral")
    for alpha in (0.05, 0.5, 1.0, 2.0, 4.0):
        x, Iu, _ = line(alpha, PhaseDrive(S=S, alpha_ac=AAC, theta_fn=uniform(Q0)), 1.5037, kin=False)
        pu = np.trapezoid(Iu, x)
        for beta in (0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4):
            d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=spatial_chirp(Q0, beta))
            x, I, K = line(alpha, d, 1.5037, qmax=Q0 * (1 + beta / 2))
            a, kap = np.pi * alpha * J1, Q0 * beta
            G = a * a / kap if beta else np.inf
            P = 1 - np.exp(-2 * np.pi * G) if beta else np.nan
            mid = np.abs(x) < 0.3 * Q0 * beta / 2 if beta else np.abs(x) < 0.5
            f = features(x, I)
            print(f" {alpha:4.2f} {beta:5.2f} | {G:10.3g} | {P:6.3f} | {np.median(I[mid]):9.3f} {f['peak']:6.3f} "
                  f"{f['width']:6.2f} {f['integral']:8.3f} {f['npk']:3d} | {np.median(K[mid]):9.3f} "
                  f"{np.trapezoid(K, x):8.3f} | {pu:8.3f}")


def section_A2():
    print("== A2 retarded-time chirp: satellite trajectory while the front fills (alpha = 0.5, 2) ==")
    for alpha in (0.5, 2.0):
        for b in (0.1, 0.3):
            d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=retarded_chirp(Q0, b))
            eta = np.arange(Q0 * 0.6, Q0 * 1.6, 0.25)
            t = np.linspace(0.1, 2.0, 20) + 0.0037
            nz = required_nz(alpha, AAC * (1 + 2 * b), S=S, eta_max=eta.max())
            I = rocking_cube(alpha, d, eta, t, n_z=nz)
            rows = []
            for k in range(t.size):
                f = features(eta - Q0, I[k])
                pred = Q0 * b * (t[k] - min(t[k], 1.0) / 2)      # mean local q over the filled depth
                rows.append(f"{t[k]:.2f}:{f['centroid']:+.1f}({pred:+.1f}),w{f['width']:.1f},p{f['npk']}")
            print(f"  alpha={alpha} b={b}: t:centroid(pred), width, npk ->\n    " + "  ".join(rows))


def section_E():
    print("== E focusing search: theta = q0 (zeta - t) + c2 (zeta-1/2)^2 + c3 (zeta-1/2)^3, t = 1.5 ==")
    for alpha in (0.05, 0.5, 1.0, 2.0):
        x, Iu, Ku = line(alpha, PhaseDrive(S=S, alpha_ac=AAC, theta_fn=uniform(Q0)), 1.5037, half=0.12, step=0.05)
        fu = features(x, Iu)
        best = []
        for c2 in np.linspace(-60, 60, 13):
            for c3 in np.linspace(-120, 120, 9):
                th = (lambda c2, c3: lambda z, t: Q0 * (z - t) + c2 * (z - 0.5) ** 2 + c3 * (z - 0.5) ** 3)(c2, c3)
                d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=th)
                x, I, K = line(alpha, d, 1.5037, half=0.12, step=0.05, kin=False, qmax=Q0 + 200)
                f = features(x, I)
                best.append((f["peak"] / f["integral"], f["width"], c2, c3, f["peak"], f["npk"]))
        r_u = fu["peak"] / fu["integral"]
        top = max(best)
        narrow = min(best, key=lambda b: b[1])
        print(f"  alpha={alpha}: uniform peak/int {r_u:.4f}, FWHM {fu['width']:.2f} | best peak/int "
              f"{top[0]:.4f} at c2={top[2]:.0f}, c3={top[3]:.0f} (ratio {top[0]/r_u:.3f}) | narrowest FWHM "
              f"{narrow[1]:.2f} at c2={narrow[2]:.0f}, c3={narrow[3]:.0f}, peaks {narrow[5]}, peak {narrow[4]:.3f}")


if __name__ == "__main__":
    for s in (sys.argv[1:] or ["A1", "A2", "E"]):
        globals()[f"section_{s}"]()
