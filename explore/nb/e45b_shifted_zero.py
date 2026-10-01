"""Entry 45b: renormalised Bessel zero. Near S = j_{n,1} the effective coupling of channel n is
g_eff = pi alpha [J_{-n}(S) - (pi alpha / q)^2 C_n(S)] (Entry 45). It is real, so tuning S alone
can make it vanish exactly at
    S* = j + (R / 2)^2 C_n(j) / J_{-n}'(j),   R = alpha / alpha_ac      (first order in R^2)
Test: scan S around j, find the minimum of the channel population max_eta |b_{-n}(1)|^2,
compare S* (numerical) with the prediction, and check how deep the minimum is.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402
from scipy.special import jv, jvp  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from e45_zero import ZEROS, C_n, channel_population  # noqa: E402

if __name__ == "__main__":
    print(" n  alpha  aac    R   | S* numeric   S* pred    j      | min population  (at j: pop)")
    for n in (0, 1, 2):
        j = ZEROS[n]
        Cn = C_n(n, j)
        dJ = jvp(-n, j) if n else jvp(0, j)
        for alpha, aac in ((1.0, 5.0), (1.0, 3.0), (2.0, 5.0), (2.0, 3.0), (1.0, 2.0), (2.0, 2.0)):
            R = alpha / aac
            f = lambda S: np.log(channel_population(alpha, aac, S, n)[0] + 1e-300)  # noqa: E731
            pred = j + (R / 2) ** 2 * Cn / dJ
            span = max(4 * abs(pred - j), 0.02)
            res = minimize_scalar(f, bounds=(j - span, j + span), method="bounded", options={"xatol": 1e-6})
            print(f" {n:+d} {alpha:5.2f} {aac:4.1f} {R:5.2f} | {res.x:.6f}  {pred:.6f}  {j:.6f} | "
                  f"{np.exp(res.fun):.2e}  ({channel_population(alpha, aac, j, n)[0]:.2e})", flush=True)
