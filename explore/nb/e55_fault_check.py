"""Entry 55: compare the Entry 39 phase-wall results with the plane-wave stacking-fault formula of
Authier & Simon, Acta Cryst. A 24, 517 (1968), eq. (4-5), non-absorbing symmetric Laue case.

Mapping: satellite n is a two-beam crystal with coupling a = pi alpha |J_n(S)| (thickness 1) and a
fault phase delta = n dphi at depth zeta_b. Authier-Simon variables: B l -> a, B eta_AS -> eta / 2,
u l -> s = sqrt(a^2 + eta^2/4), l_I = zeta_b, l_II = 1 - zeta_b (which side is first does not change the
intensity). Eq. (4-5), with sin(lu)/u normalised so that delta = 0 gives the perfect-crystal amplitude:

  A = a e^{i delta/2} { cos(d/2) sin s / s - i sin(d/2) sin(s (l_I - l_II)) / s
                        + (eta/2) / s^2 sin(d/2) [cos(s (l_I - l_II)) - cos s] }

Checked here against an independent two-segment transfer-matrix product (eq. 4-3 structure).
Usage: python e55_fault_check.py
"""
import numpy as np
from scipy.linalg import expm
from scipy.special import jv

AAC, S = 40.0, 1.5
Q = 2 * np.pi * AAC
ETA = np.arange(-0.45 * Q, 0.45 * Q, 0.1)          # Entry 39 window around n q, step 0.1


def as_eq45(a, eta, delta, zb):
    s = np.sqrt(a * a + eta * eta / 4)
    dl = zb - (1 - zb)
    return a * np.exp(1j * delta / 2) * (np.cos(delta / 2) * np.sin(s) / s - 1j * np.sin(delta / 2) * np.sin(s * dl) / s
                                         + (eta / 2) / s**2 * np.sin(delta / 2) * (np.cos(s * dl) - np.cos(s)))


def matrix(a, eta, delta, zb):
    out = np.empty(eta.size, complex)
    for i, e in enumerate(eta):
        H1 = np.array([[-e / 2, a], [a, e / 2]])
        H2 = np.array([[-e / 2, a * np.exp(1j * delta)], [a * np.exp(-1j * delta), e / 2]])
        out[i] = (expm(-1j * (1 - zb) * H2) @ expm(-1j * zb * H1) @ np.array([1, 0]))[1]
    return out


def integral(a, delta, zb):
    return np.trapezoid(np.abs(as_eq45(a, ETA, delta, zb)) ** 2, ETA)


if __name__ == "__main__":
    e = np.linspace(-30, 30, 13)
    for a, d, zb in ((3.5, np.pi, 0.5), (3.5, np.pi / 2, 0.3), (1.46, 2.0, 0.7)):
        dev = np.max(np.abs(np.abs(as_eq45(a, e, d, zb)) ** 2 - np.abs(matrix(a, e, d, zb)) ** 2))
        print(f"eq. (4-5) vs transfer matrix: a={a}, delta={d:.3f}, zeta_b={zb}: max |dI| = {dev:.1e}")

    print("\nEntry 39 C (wall at 0.5): integral(dphi) / integral(0), TT (Entry 39) vs plane-wave fault formula")
    tt = {(0.5, 1): [3.70, 3.77, 3.96, 4.14, 4.19, 3.94], (2.0, 1): [10.1, 10.8, 12.1, 13.0, 13.4, 11.8],
          (2.0, 2): [6.49, 7.88, 9.45, 8.09, 6.50, 9.43], (2.0, 3): [0.85, 0.88, 0.88, 0.85, 0.93, 0.98]}
    dphis = np.array([0, 0.25, 0.5, 0.75, 1.0, 1.5]) * np.pi
    for (alpha, n), v in tt.items():
        a = np.pi * alpha * abs(jv(n, S))
        f = [integral(a, n * dp, 0.5) for dp in dphis]
        print(f" alpha={alpha} n={n} a={a:.3f}: formula abs {f[0]:.3f} | ratios TT " +
              " ".join(f"{x / v[0]:.3f}" for x in v) + " | formula " + " ".join(f"{x / f[0]:.3f}" for x in f))

    print("\nEntry 39 Cz (dphi = pi, n = 1): ratio pi-wall / no wall vs zeta_b")
    for alpha in (0.5, 2.0):
        a = np.pi * alpha * abs(jv(1, S))
        print(f" alpha={alpha}: " + " ".join(f"{zb:.1f}:{integral(a, np.pi, zb) / integral(a, 0.0, zb):.2f}"
                                             for zb in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)))
        e0 = np.array([0.0])
        print("   centre abs(A)/a vs a abs(1-2 zeta_b): " + " ".join(
            f"{zb:.2f}:{abs(as_eq45(a, e0, np.pi, zb)[0]):.3f}/{abs(np.sin(a * (1 - 2 * zb))):.3f}"
            for zb in (0.05, 0.28, 0.5, 0.72)))
