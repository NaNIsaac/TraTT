"""Entry 40: coarse automated anomaly search over simple acoustic profiles.

Families (n = 1, snapshot t = 1.5, fully filled): spatial chirp beta; q-wall dc at depth zb;
phase wall dphi at depth zb; amplitude ramp T (non-uniform amplitude). alpha in
{0.05, 0.5, 1, 2, 4}, S in {1.0, 2.5}. For every case the full-TT satellite line (static
subtracted) and the kinematic line are reduced to features and compared with the
uniform-q control at the same (alpha, S). Scores:
  topo   |npk - npk_control|                         (peak creation/destruction)
  dint   |ln(integral / integral_control)|           (spectral weight created/destroyed)
  dyn    dint(TT) - dint(kinematic)                  (part not explained by Fourier optics)
  narrow ln(width_control / width) if narrower       (narrow features)
  conc   ln[(peak/int) / (peak/int)_control]         (concentration of weight)
Output: ranked table (top 15 by dyn, then by topo+narrow+conc), saved to e40_search.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from qmod import PhaseDrive, features, kinematic, phase_wall, q_wall, spatial_chirp, uniform  # noqa: E402
from tratt.core import required_nz, rocking  # noqa: E402

AAC = 40.0
Q0 = 2 * np.pi * AAC
T_SNAP = 1.5037


def profiles():
    for beta in (0.02, 0.05, 0.1, 0.2, 0.4):
        yield ("chirp", f"beta={beta}", dict(theta_fn=spatial_chirp(Q0, beta)), Q0 * (1 + beta / 2))
    for dc in (0.3, 1.0, 2.0, 4.0):
        for zb in (0.25, 0.5, 0.75):
            yield ("qwall", f"dc={dc},zb={zb}", dict(theta_fn=q_wall(Q0, Q0 + 2 * np.pi * dc, T_SNAP - zb)),
                   Q0 + 2 * np.pi * dc)
    for dphi in (np.pi / 2, np.pi):
        for zb in (0.25, 0.5, 0.75):
            yield ("pwall", f"dphi={dphi/np.pi:.2f}pi,zb={zb}",
                   dict(theta_fn=phase_wall(Q0, dphi, T_SNAP - zb)), Q0)
    for T in (0.3, 1.0):
        yield ("ramp", f"T={T}", dict(theta_fn=uniform(Q0), amp=lambda s, T=T: np.clip(s / T, 0, 1)), Q0)


def line(alpha, S, kw, qmax, half=0.45, step=0.1, nzmul=1):
    eta = np.arange(Q0 * (1 - half), Q0 * (1 + half) + (qmax - Q0), step)
    nz = nzmul * required_nz(alpha, qmax / (2 * np.pi), S=S, eta_max=eta.max())
    d = PhaseDrive(S=S, alpha_ac=AAC, **kw)
    I = rocking(alpha, eta, phi=d.phi(T_SNAP), n_z=nz) - rocking(alpha, eta, n_z=nz)
    K = kinematic(alpha, d, eta, T_SNAP) - kinematic(alpha, PhaseDrive(S=0.0, alpha_ac=AAC), eta, T_SNAP)
    return eta - Q0, I, K


def score(f, c):
    narrow = max(0.0, np.log(c["width"] / f["width"])) if f["width"] > 0 else 0.0
    return dict(topo=abs(f["npk"] - c["npk"]), dint=abs(np.log(f["integral"] / c["integral"])),
                narrow=narrow, conc=np.log((f["peak"] / f["integral"]) / (c["peak"] / c["integral"])))


def main():
    rows = []
    for alpha in (0.05, 0.5, 1.0, 2.0, 4.0):
        for S in (1.0, 2.5):
            x, Ic, Kc = line(alpha, S, dict(theta_fn=uniform(Q0)), Q0)
            fc, fkc = features(x, Ic), features(x, Kc)
            for fam, lab, kw, qmax in profiles():
                x, I, K = line(alpha, S, kw, qmax)
                f, fk = features(x, I), features(x, K)
                s, sk = score(f, fc), score(fk, fkc)
                rows.append(dict(fam=fam, lab=lab, alpha=alpha, S=S, **{k: float(v) for k, v in s.items()},
                                 dyn=float(s["dint"] - sk["dint"]), kin_dint=float(sk["dint"]),
                                 int_ratio=float(f["integral"] / fc["integral"]), npk=f["npk"], npk_c=fc["npk"],
                                 width=f["width"], width_c=fc["width"], peak=f["peak"], peak_c=fc["peak"]))
            print(f"done alpha={alpha} S={S}", flush=True)
    json.dump(rows, open(HERE / "e40_search.json", "w"), indent=1)
    hdr = " fam   label                  alpha  S   | int/ctrl  dyn   kin_dint | npk/ctrl | width/ctrl | peak/ctrl | conc"
    for key, title in (("dyn", "by dynamical weight change"),
                       ("shape", "by topology + narrowing + concentration")):
        print(f"== top 15 {title} ==\n{hdr}")
        k = (lambda r: r["dyn"]) if key == "dyn" else (lambda r: r["topo"] + 2 * r["narrow"] + max(r["conc"], 0))
        for r in sorted(rows, key=k, reverse=True)[:15]:
            print(f" {r['fam']:5s} {r['lab']:22s} {r['alpha']:4.2f} {r['S']:3.1f} | {r['int_ratio']:7.3f} "
                  f"{r['dyn']:6.3f} {r['kin_dint']:6.3f} | {r['npk']:2d}/{r['npk_c']:2d} | "
                  f"{r['width']:5.1f}/{r['width_c']:5.1f} | {r['peak']:.3f}/{r['peak_c']:.3f} | {r['conc']:+.2f}")


if __name__ == "__main__":
    main()
