"""Entry 65 follow-up: is a few-cycle *packet* linear where Entry 65 says it should be?

Entry 65 estimated the packet's satellite amplitude as nu ~ pi alpha |J_1(S)| tau_d (tau_d = packet
decay length in units of L, v = 1) and the first-order error as ~nu^2. Test it with full TT on the
Entry 59 packet (e58_chalise_ref.run): alpha = 40, kappa = 6 (g = 0.024), R_D = 0.2, and
S = 0.02, 0.2, 0.7. The weight W / S against the S = 0.02 result measures the nonlinearity.
Usage: python e65_packet.py
"""
import json

import numpy as np
from scipy.special import jv

from e58_chalise_ref import packet, run

alpha, kappa, RD = 40.0, 6.0, 0.2
wd = 2 * np.pi * alpha
q0 = wd / RD
eta = q0 * np.array([0.5, 0.75, 1.0, 1.25, 1.5, 2.0])
_, tau_d = packet(1.0, q0)
out = {}
base = None
for S in (0.02, 0.2, 0.7):
    r = run(alpha, S, q0, eta, kappa=kappa)
    w = np.asarray(r["W"]) / S
    base = w if base is None else base
    nu = np.pi * alpha * abs(jv(1, S)) * tau_d
    dev = np.max(np.abs(w / base - 1))
    print(f"S={S:<5} nu_est={nu:.3f} nu^2={nu**2:.3f}  max |W/S / (W/S)_lin - 1| = {dev:.3f}   "
          f"per eta: " + " ".join(f"{v:.3f}" for v in w / base), flush=True)
    out[str(S)] = dict(nu=nu, dev=float(dev), ratio=(w / base).tolist(), n_z=int(r["nz"]))
json.dump(dict(alpha=alpha, kappa=kappa, RD=RD, eta=eta.tolist(), tau_d=tau_d, runs=out),
          open("e65_packet.json", "w"), indent=1)
