"""Entry 30: history-dependent diffraction.

(a) Triangular pulse 0 -> S_max -> 0 over 2T at the transducer. Compare I_n when the
    transducer amplitude S A(t) is the same on the way up and on the way down.
(b) After the transducer is silent (t > 2T) the packet is still inside the crystal.
(c) Multiple revival: slow triangle to S_max = 10 sweeps J_1 through its zeros twice.
(d) Phase-coded packet: box of length 0.8 with a pi carrier flip at 0.4. While both
    halves are inside, odd orders cancel and even orders add.
Laue and Bragg, alpha = 0.5 (no Pendellosung). Theory: Theta_n read by sin^2 / tanh^2.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from memory_drive import HistoryDrive, box, theta, triangle  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

READ = {"laue": lambda x: np.sin(x) ** 2, "bragg": lambda x: np.tanh(x) ** 2}
ALPHA, AAC = 0.5, 40.0


def run(d, t, ns, geom):
    I = rocking_cube(ALPHA, d, np.array([n * d.q for n in ns]), t, geometry=geom, n_z=2400)
    th = np.stack([READ[geom](np.pi * ALPHA * np.abs(theta(d, n, t))) for n in ns], axis=1)
    return I, th


print("(a)+(b) triangular pulse, S_max = 3")
for geom in ("laue", "bragg"):
    for T in (0.25, 1.0, 4.0):
        d = HistoryDrive(S=3.0, alpha_ac=AAC, env=triangle(T))
        a = np.array([0.3, 0.6, 0.9])                        # A on the way up and down
        tu, td = a * T, (2 - a) * T
        I, th = run(d, np.concatenate([tu, td, [2 * T + 0.5 * min(1.0, 2 * T)]]) + 0.0037, (1, 2), geom)
        rows = " ".join(f"A={x:.1f}: n1 {I[i,0]:.3f}/{I[3+i,0]:.3f} n2 {I[i,1]:.3f}/{I[3+i,1]:.3f};"
                        for i, x in enumerate(a))
        print(f"  {geom:5s} T={T}: up/down {rows}  after drive (A=0): n1 {I[6,0]:.3f} n2 {I[6,1]:.3f}"
              f" | law err {np.max(np.abs(I-th)):.3f}")

print("(c) multiple revival, S_max = 10, T = 3 (slow) and T = 0.5 (fast), n = 1")
for T in (3.0, 0.5):
    d = HistoryDrive(S=10.0, alpha_ac=AAC, env=triangle(T))
    t = np.linspace(0.0, 2 * T + 1.0, 241) + 0.0037
    I, th = run(d, t, (1,), "laue")
    y = I[:, 0]
    pk = [i for i in range(1, y.size - 1) if y[i] >= y[i - 1] and y[i] > y[i + 1] and y[i] > 0.05]
    dk = [i for i in range(1, y.size - 1) if y[i] <= y[i - 1] and y[i] < y[i + 1] and y[i] < 0.02]
    print(f"  T={T}: bright maxima at t={np.round(t[pk], 2)} values {np.round(y[pk], 3)}; "
          f"dark minima at t={np.round(t[dk], 2)}; law err {np.max(np.abs(y - th[:, 0])):.3f}")

print("(d) phase-coded packet: box 0.8, pi flip at 0.4, S = 1.5")
flip = lambda s: np.where(s >= 0.4, 1.0, 0.0)  # noqa: E731
for geom in ("laue", "bragg"):
    d = HistoryDrive(S=1.5, alpha_ac=AAC, env=box(0.8), flip=flip)
    t = np.array([0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.9]) + 0.0037
    I, th = run(d, t, (1, 2), geom)
    print(f"  {geom:5s} t:  " + " ".join(f"{x:5.2f}" for x in t))
    print(f"        n=1: " + " ".join(f"{v:5.3f}" for v in I[:, 0]))
    print(f"        n=2: " + " ".join(f"{v:5.3f}" for v in I[:, 1]) + f"   law err {np.max(np.abs(I-th)):.3f}")
