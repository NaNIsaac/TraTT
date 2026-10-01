"""Entry 39: acoustic domain walls — q-jump (B) and phase jump (C).

Transducer emits q1 (or phase 0) for one transit, then q2 (or phase dphi). For t in [1, 2]
the wall sits at depth zeta_b = t - 1: shallow side new (q2 / dphi), deep side old.
Controls for every case:
  incoherent  (I_deep_alone - I0) + (I_shallow_alone - I0), the other side unmodulated
  kinematic   first-Born Fourier model of the same profile
  smoothed    the wall smoothed over one acoustic wavelength (w = 1 / alpha_ac in s)
Usage: python e39_walls.py [B] [Bt] [C] [Cz]
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from qmod import PhaseDrive, features, kinematic, phase_wall, q_wall  # noqa: E402
from tratt.core import required_nz, rocking  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

AAC, S = 40.0, 1.5
Q1 = 2 * np.pi * AAC
W1 = 1.0 / AAC  # one acoustic period in retarded time = one wavelength in depth


def amp_old(s):
    return np.where((s >= 0) & (s < 1.0), 1.0, 0.0)


def amp_new(s):
    return np.where(s >= 1.0, 1.0, 0.0)


def sat(alpha, theta_fn, eta, t, amp=None, nz=None):
    d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=theta_fn, amp=amp)
    nz = nz or required_nz(alpha, AAC * 1.3, S=S, eta_max=eta.max())
    I0 = rocking(alpha, eta, n_z=nz)
    return rocking(alpha, eta, phi=d.phi(t), n_z=nz) - I0, d


def coherence(I, Iinc, x):
    return np.trapezoid(np.abs(I - Iinc), x) / np.trapezoid(np.abs(Iinc), x)


def section_B():
    print("== B q-wall at depth 0.5 (t = 1.5): q2 = q1 + 2 pi dc, n = 1 ==")
    print(" alpha   dc | npk both/inc/kin | integral both/inc | coh. dev |TT-inc|/inc | "
          "mid-point I both/inc | TT vs kin shape dev | smoothed: npk, integral")
    t = 1.5037
    for alpha in (0.05, 0.5, 2.0):
        for dc in (0.1, 0.3, 1.0, 2.0, 4.0, 10.0):
            q2 = Q1 + 2 * np.pi * dc
            eta = np.arange(Q1 - 60, q2 + 60, 0.1)
            x = eta - Q1
            th = q_wall(Q1, q2, 1.0)
            Ib, d = sat(alpha, th, eta, t)
            Id, _ = sat(alpha, th, eta, t, amp=amp_old)
            In, _ = sat(alpha, th, eta, t, amp=amp_new)
            Iinc = Id + In
            K = kinematic(alpha, d, eta, t) - kinematic(alpha, PhaseDrive(S=0.0, alpha_ac=AAC), eta, t)
            Is, _ = sat(alpha, q_wall(Q1, q2, 1.0, W1), eta, t)
            fb, fi, fk, fs = (features(x, y) for y in (Ib, Iinc, K, Is))
            km = np.argmin(np.abs(x - np.pi * dc))
            kdev = np.trapezoid(np.abs(Ib / fb["integral"] - K / fk["integral"]), x)
            print(f" {alpha:4.2f} {dc:5.1f} | {fb['npk']:2d}/{fi['npk']:2d}/{fk['npk']:2d} | "
                  f"{fb['integral']:7.3f}/{fi['integral']:7.3f} | {coherence(Ib, Iinc, x):6.3f} | "
                  f"{Ib[km]:.4f}/{Iinc[km]:.4f} | {kdev:6.3f} | {fs['npk']:2d}, {fs['integral']:7.3f}")


def section_Bt():
    print("== B moving wall: dc = 1 (comparable to the width), alpha = 0.5 and 2; features vs wall depth ==")
    for alpha in (0.5, 2.0):
        dc = 1.0
        q2 = Q1 + 2 * np.pi * dc
        eta = np.arange(Q1 - 60, q2 + 60, 0.1)
        x = eta - Q1
        d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=q_wall(Q1, q2, 1.0))
        t = 1.0 + np.linspace(0.0, 1.0, 21) + 0.0037
        nz = required_nz(alpha, AAC * 1.3, S=S, eta_max=eta.max())
        I = rocking_cube(alpha, d, eta, t, n_z=nz) - rocking(alpha, eta, n_z=nz)[None, :]
        print(f"  alpha={alpha}: zeta_b: npk, centroid, width, peak, integral")
        print("   " + "  ".join(f"{t[k]-1:.2f}:{f['npk']},{f['centroid']:+.1f},{f['width']:.1f},"
                               f"{f['peak']:.3f},{f['integral']:.2f}"
                               for k, f in enumerate(features(x, I[k]) for k in range(t.size))))


def section_C():
    print("== C phase wall at depth 0.5 (t = 1.5), same q: integrated satellite per order vs dphi ==")
    t = 1.5037
    for alpha in (0.05, 0.5, 2.0):
        for n in (1, 2, 3):
            eta = np.arange(n * Q1 - 0.45 * Q1, n * Q1 + 0.45 * Q1, 0.1)
            x = eta - n * Q1
            rows = []
            for dphi in (0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi, 3 * np.pi / 2):
                I, d = sat(alpha, phase_wall(Q1, dphi, 1.0), eta, t)
                K = kinematic(alpha, d, eta, t) - kinematic(alpha, PhaseDrive(S=0.0, alpha_ac=AAC), eta, t)
                f, fk = features(x, I), features(x, K)
                rows.append(f"{dphi/np.pi:.2f}pi: int {f['integral']:.3f} (kin {fk['integral']:.3f}) "
                            f"ctr {f['centre']:.3f} npk {f['npk']}")
            print(f"  alpha={alpha} n={n}: " + " | ".join(rows))
    print("  smoothing check (alpha 2, n 1, dphi pi): sharp vs w = 1 wavelength")
    eta = np.arange(0.55 * Q1, 1.45 * Q1, 0.1)
    for w in (0.0, W1, 3 * W1):
        I, _ = sat(2.0, phase_wall(Q1, np.pi, 1.0, w), eta, t)
        f = features(eta - Q1, I)
        print(f"   w={w:.4f}: integral {f['integral']:.3f}, centre {f['centre']:.4f}, npk {f['npk']}")


def section_Cz():
    print("== C stacking-fault analogue: centre and integral vs wall depth, dphi = pi, n = 1 ==")
    for alpha in (0.5, 2.0):
        eta = np.arange(0.55 * Q1, 1.45 * Q1, 0.1)
        d = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=phase_wall(Q1, np.pi, 1.0))
        d0 = PhaseDrive(S=S, alpha_ac=AAC, theta_fn=phase_wall(Q1, 0.0, 1.0))
        t = 1.0 + np.linspace(0.0, 1.0, 21) + 0.0037
        nz = required_nz(alpha, AAC, S=S, eta_max=eta.max())
        I0 = rocking(alpha, eta, n_z=nz)[None, :]
        I = rocking_cube(alpha, d, eta, t, n_z=nz) - I0
        Iu = rocking_cube(alpha, d0, eta, t, n_z=nz) - I0
        k0 = np.argmin(np.abs(eta - Q1))
        p, pu = np.trapezoid(I, eta, axis=1), np.trapezoid(Iu, eta, axis=1)
        print(f"  alpha={alpha}: zeta_b: centre(pi wall)/centre(no wall), integral ratio pi/none")
        print("   " + "  ".join(f"{t[k]-1:.2f}:{I[k,k0]:.3f}/{Iu[k,k0]:.3f},{p[k]/pu[k]:.3f}" for k in range(t.size)))


if __name__ == "__main__":
    for s in (sys.argv[1:] or ["B", "Bt", "C", "Cz"]):
        globals()[f"section_{s}"]()
