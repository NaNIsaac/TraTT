"""Entry 31: does the brightest satellite order switch in time?

Linear ramp to S = 6 over T = 4, alpha = 0.5 (weak), Laue, alpha_ac = 40. Track
argmax_n I(nq, t) over n = 1..5; compare with (i) the history law and (ii) the
instantaneous (quasi-static) Bessel ladder J_n(S A(t)), to measure the memory lag.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from memory_drive import HistoryDrive, ramp, theta  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

ALPHA, NS = 0.5, (1, 2, 3, 4, 5)
for T in (4.0, 1.0):
    d = HistoryDrive(S=6.0, alpha_ac=40.0, env=ramp(T))
    t = np.linspace(0.0, T + 1.5, 301) + 0.0037
    I = rocking_cube(ALPHA, d, np.array([n * d.q for n in NS]), t)  # n_z from the extended rule
    law = np.stack([np.sin(np.pi * ALPHA * np.abs(theta(d, n, t))) ** 2 for n in NS], 1)
    inst = np.stack([np.sin(np.pi * ALPHA * np.abs(jv(n, 6.0 * d.envelope(t)))) ** 2 for n in NS], 1)

    def switches(M):
        lead = np.array(NS)[np.argmax(M, 1)]
        lead[M.max(1) < 0.02] = 0
        ch = [(round(float(t[i]), 2), int(lead[i - 1]), int(lead[i])) for i in range(1, t.size)
              if lead[i] != lead[i - 1]]
        return ch

    print(f"T={T}: law err {np.max(np.abs(I - law)):.3f}")
    print("  brightest-order switches (t, from, to), full TT:", switches(I))
    print("  history law                               :", switches(law))
    print("  instantaneous J_n(S A(t)) (no memory)     :", switches(inst))
