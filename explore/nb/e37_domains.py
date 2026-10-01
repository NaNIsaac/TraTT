"""Entry 37: two acoustic domains that each diffract strongly but cancel together.

Transducer: amplitude S2 for 0 <= s < 1, then S1. For t in [1, 2] the crystal holds S1 in
zeta < z_b = t - 1 (entrance side, recent) and S2 deeper; the domain wall moves at v_s.
Law: Theta_n = z_b J_n(S1) + (1 - z_b) J_n(S2), zero at z_b = r* = J_n(S2) / (J_n(S2) - J_n(S1)).
Optimum for n = 1: S1 = 1.841 (J_1 max), S2 = 5.331 (J_1 min) -> r* = 0.373.
Checks: each domain alone (the other region unmodulated), full line shape and window
integral at the dark moment, other orders (order selection), Bragg, sharpness of the dark.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from memory_drive import HistoryDrive, theta  # noqa: E402
from tratt.core import required_nz, rocking  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

AAC = 40.0
S1 = minimize_scalar(lambda s: -jv(1, s), bounds=(1, 3), method="bounded").x
S2 = minimize_scalar(lambda s: jv(1, s), bounds=(4, 7), method="bounded").x
R = jv(1, S2) / (jv(1, S2) - jv(1, S1))
print(f"S1={S1:.4f} (J1={jv(1,S1):+.4f}), S2={S2:.4f} (J1={jv(1,S2):+.4f}), r*={R:.4f}")


def drive(a2, a1):
    return HistoryDrive(S=1.0, alpha_ac=AAC, env=lambda s: np.where(s < 1.0, a2, a1))


READ = {"laue": lambda x: np.sin(x) ** 2, "bragg": lambda x: np.tanh(x) ** 2}
for geom in ("laue", "bragg"):
    for alpha in (0.5, 2.0):
        nz = required_nz(alpha, AAC, S=S2, eta_max=3 * 2 * np.pi * AAC)
        d = drive(S2, S1)
        t = np.linspace(0.0, 2.3, 231) + 0.0037
        I = rocking_cube(alpha, d, np.array([n * d.q for n in (1, 2, 3)]), t, geometry=geom, n_z=nz)
        law = np.stack([READ[geom](np.pi * alpha * np.abs(theta(d, n, t))) for n in (1, 2, 3)], 1)
        ts = 1.0 + R + 0.0037
        both = rocking(alpha, np.array([d.q]), phi=d.phi(ts), geometry=geom, n_z=nz)[0]
        a1 = rocking(alpha, np.array([d.q]), phi=drive(0.0, S1).phi(ts), geometry=geom, n_z=nz)[0]
        a2 = rocking(alpha, np.array([d.q]), phi=drive(S2, 0.0).phi(ts), geometry=geom, n_z=nz)[0]
        m = (t > 1.0) & (t < 2.0)
        k = np.argmin(np.where(m, I[:, 0], 9))
        ok = m & (I[:, 0] < 0.01 * I[m, 0].max())
        width = (t[ok][-1] - t[ok][0]) if ok.any() else 0.0
        print(f"--- {geom} alpha={alpha}: law max err (n=1,2,3) {np.round(np.abs(I-law).max(0), 4)}")
        print(f"    at the wall position r*: S1 domain alone {a1:.3f}, S2 domain alone {a2:.3f}, "
              f"both {both:.2e}")
        print(f"    n=1 minimum {I[k,0]:.2e} at t={t[k]:.3f} (pred {1+R:.3f}); <1% of max for dt={width:.3f} t_tr; "
              f"n=2, n=3 there: {I[k,1]:.3f}, {I[k,2]:.3f}")

# full line shape at the dark moment vs single-domain dark (Pendellosung zero) for comparison
alpha = 2.0
d = drive(S2, S1)
q = d.q
eta = np.linspace(0.5 * q, 1.5 * q, 2001)
nz = required_nz(alpha, AAC, S=S2, eta_max=1.5 * q)
I_dark = rocking(alpha, eta, phi=d.phi(1.0 + R + 0.0037), n_z=nz)
I_S1 = rocking(alpha, eta, phi=drive(0.0, S1).phi(1.0 + R + 0.0037), n_z=nz)
p = lambda y: np.trapezoid(y, eta)  # noqa: E731
k0 = np.argmin(np.abs(eta - q))
pk = [round(float(eta[i] - q), 2) for i in range(1, eta.size - 1)
      if I_dark[i] > I_dark[i - 1] and I_dark[i] >= I_dark[i + 1] and I_dark[i] > 0.05 * I_dark.max()]
print(f"line shape at the cancellation (alpha=2): centre {I_dark[k0]:.2e}, window integral {p(I_dark):.3f} "
      f"(S1 domain alone: centre {I_S1[k0]:.3f}, integral {p(I_S1):.3f}); max {I_dark.max():.3f}; "
      f"side maxima at eta-q = {pk[:8]}")

# order selection: cancel n = 1 while keeping n = 2 (or the reverse), weak coupling alpha = 0.5
print("order selection (domain fractions at which each order cancels):")
for Sa, Sb in ((S1, S2), (1.0, 4.5), (2.5, 4.2), (3.0, 6.5)):
    rs = {}
    for n in (1, 2, 3):
        ja, jb = jv(n, Sa), jv(n, Sb)
        rs[n] = jb / (jb - ja) if ja * jb < 0 else None
    print(f"  S_a={Sa:.3f} S_b={Sb:.3f}: r* for n=1,2,3 = " +
          ", ".join(f"{rs[n]:.3f}" if rs[n] is not None else "none" for n in (1, 2, 3)) +
          " | Theta_2, Theta_3 at r*_1: " +
          (", ".join(f"{rs[1]*jv(n,Sa)+(1-rs[1])*jv(n,Sb):+.3f}" for n in (2, 3)) if rs[1] else "-"))
