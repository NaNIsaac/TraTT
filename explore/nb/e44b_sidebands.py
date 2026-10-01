"""Entry 44b: how many forward-beam Floquet sidebands a_k must be kept, vs R = alpha / alpha_ac?

All b_m are kept; a_k is truncated to |k| <= kmax. k_req = smallest kmax with E_pop < 0.02
against the converged solution. Also the analogous truncation of the diffracted harmonics
(all a_k kept, |m + n_dom| limited) is reported for comparison. Output feeds the
convergence figure.
"""
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from tratt.floquet import FloquetLaue, auto_K  # noqa: E402
from e44_mechanism import windows  # noqa: E402

if __name__ == "__main__":
    out = []
    for S in (1.5, 4.0):
        for aac in (5.0, 2.0):
            for R in (0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
                alpha = R * aac
                q = 2 * np.pi * aac
                N = int(np.ceil(S)) + 4
                half = N * q + 4 * np.pi * alpha + 20
                eta = np.arange(-half, half, 0.25)
                K = max(auto_K(S), N) + 3 + int(np.ceil(4 * np.pi * alpha / q))
                ns = np.arange(-N, N + 1)
                pf = windows(eta, FloquetLaue(alpha, aac, S, K).rocking(eta), q, ns)
                errs = []
                for kmax in range(0, 9):
                    I = FloquetLaue(alpha, aac, S, K, mask=lambda k, m, kmax=kmax: abs(k) <= kmax).rocking(eta)
                    errs.append(float(np.abs(windows(eta, I, q, ns) - pf).sum() / pf.sum()))
                kreq = next((i for i, e in enumerate(errs) if e < 0.02), None)
                out.append(dict(S=S, aac=aac, R=R, errs=errs, kreq=kreq))
                print(f"S={S} aac={aac} R={R}: k_req={kreq}  errs " + " ".join(f"{e:.3f}" for e in errs), flush=True)
    json.dump(out, open(HERE / "e44b_sidebands.json", "w"))
