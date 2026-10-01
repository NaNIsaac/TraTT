"""Entry 28: transient-only satellite. S_f at the first zero of J_1, linear ramp of length T.

I_1 = 0 before the front, I_1 = 0 once the ramp has left the crystal (Theta_1 = J_1(S_f) = 0),
I_1 > 0 in between. Laue and Bragg, alpha = 0.5 (no Pendellosung) and 2.
Also checks item 3: is n = 1 transiently the strongest order although it is absent
in the final pattern?
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jn_zeros, jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from memory_drive import HistoryDrive, ramp, theta  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

SF = jn_zeros(1, 1)[0]
READ = {"laue": lambda x: np.sin(x) ** 2, "bragg": lambda x: np.tanh(x) ** 2}

for geom in ("laue", "bragg"):
    for alpha in (0.5, 2.0):
        for T in (0.25, 0.5, 1.0, 2.0, 4.0):
            d = HistoryDrive(S=SF, alpha_ac=40.0, env=ramp(T))
            t = np.linspace(0.0, 1.0 + T + 0.3, 161) + 0.0037
            I = rocking_cube(alpha, d, np.array([d.q, 2 * d.q, 3 * d.q]), t, geometry=geom, n_z=2400)
            th1 = READ[geom](np.pi * alpha * np.abs(theta(d, 1, t)))
            I1 = I[:, 0]
            k = int(np.argmax(I1))
            above = t[I1 > 0.5 * I1[k]]
            life = above[-1] - above[0] if above.size else 0.0
            dom = np.mean(I1 > np.maximum(I[:, 1], I[:, 2]))  # fraction of samples where n = 1 leads
            print(f"{geom:5s} alpha={alpha} T={T}: max I1={I1[k]:.3f} at t={t[k]:.2f} (law {th1.max():.3f}) "
                  f"| FWHM life={life:.2f} | final I1={I1[-1]:.1e}, I2={I[-1,1]:.3f}, I3={I[-1,2]:.3f} "
                  f"| law err={np.max(np.abs(I1-th1)):.3f} | n=1 strongest in {100*dom:.0f}% of samples")
