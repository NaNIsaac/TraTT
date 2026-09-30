"""Entry 16: where does the filled-depth law break? Error of the centre law vs alpha / alpha_ac."""
import sys
from pathlib import Path

import numpy as np
from scipy.special import jv

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from tratt.cube import rocking_cube  # noqa: E402
from tratt.drive import Drive  # noqa: E402

S, n = 1.5, 1
t = np.linspace(0.0, 1.0, 101) + 0.0037
print("alpha_ac  alpha  alpha/alpha_ac  max|I-sin2|  rms")
rows = []
for aac in (5.0, 10.0, 20.0, 40.0):
    d = Drive(S=S, alpha_ac=aac)
    for alpha in (0.5, 1.0, 2.0, 4.0, 8.0):
        I = rocking_cube(alpha, d, np.array([n * d.q]), t, n_z=int(max(1500, 60 * max(aac, alpha))))[:, 0]
        e = I - np.sin(np.pi * alpha * abs(jv(n, S)) * t) ** 2
        rows.append((alpha / aac, np.max(np.abs(e))))
        print(f"{aac:8.0f} {alpha:6.1f} {alpha/aac:14.3f} {np.max(np.abs(e)):12.3f} {np.sqrt(np.mean(e**2)):6.3f}")
r = np.array(rows)
ok = r[:, 0] <= 0.1
print(f"worst max error with alpha/alpha_ac <= 0.1: {r[ok, 1].max():.3f};  with 0.1 < ratio <= 0.2: "
      f"{r[(r[:, 0] > 0.1) & (r[:, 0] <= 0.2), 1].max():.3f}")
