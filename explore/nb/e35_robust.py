"""Entry 35: robustness / failure boundary of the transient-only satellite (n = 1, S_f = j_11).

Observables (Laue centre intensity at eta = q unless stated):
  I_max   max over the transient
  I_fin   mean over one acoustic period at the end of the window (avoids ripple nodes)
  C       I_max / I_fin  (contrast)
Cases: FAST  alpha = 1,   T = 0.5  (T < t_tr: the genuine transient)
       MID   alpha = 0.5, T = 2    (Entry 28 baseline)
       STEP  alpha = 1,   T = 0    (control: ideal I_1 = 0 at all times)
       SLOW  alpha = 0.5, T = 8    (adiabatic scan; reference only)
Usage: python e35_robust.py [section ...]   sections: base atten disp nl ring combo retune integ
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import jn_zeros, jv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from nonideal import GeneralDrive, dispersed_ramp, ramp, ringing_ramp, theta_law  # noqa: E402
from tratt.core import required_nz  # noqa: E402
from tratt.cube import rocking_cube  # noqa: E402

Z = jn_zeros(1, 1)[0]
AAC = 40.0
CASES = {"FAST": (1.0, 0.5), "MID": (0.5, 2.0), "STEP": (1.0, 0.0), "SLOW": (0.5, 8.0)}


def nz_for(d, extra=0.0):
    return required_nz(1.0, d.alpha_ac, S=d.S * (1 + extra) + 2 * d.eps2 * d.S, eta_max=d.q)


def run(d, alpha, T, extra_t=0.5, nz=None, eta=None, law=True):
    """Centre (or eta-integrated if eta is an array) I_1 vs t; returns dict of observables."""
    t_end = 1.0 + T + extra_t
    t = np.concatenate([np.linspace(0.0, t_end - 1 / AAC, 81) + 0.0037,
                        np.linspace(t_end - 1 / AAC, t_end, 16)])
    nz = nz or nz_for(d)
    if eta is None:
        I = rocking_cube(alpha, d, np.array([d.q]), t, n_z=nz)[:, 0]   # all times in one sweep
    else:
        I = np.trapezoid(rocking_cube(alpha, d, eta, t, n_z=nz), eta, axis=1)
    out = dict(I_max=I[:81].max(), t_max=t[np.argmax(I[:81])], I_fin=I[81:].mean(), I=I, t=t)
    out["C"] = out["I_max"] / max(out["I_fin"], 1e-12)
    if law and eta is None:
        th = np.sin(np.pi * alpha * np.abs(theta_law(d, 1, t))) ** 2
        out["law_err"] = np.max(np.abs(I - th))
        out["law_fin"] = th[81:].mean()
    return out


def fmt(o):
    s = f"I_max {o['I_max']:.4f} (t {o['t_max']:.2f}) | I_fin {o['I_fin']:.2e} | C {o['C']:9.1f}"
    if "law_err" in o:
        s += f" | law err {o['law_err']:.4f}, law I_fin {o['law_fin']:.2e}"
    return s


def section_base():
    print("== baseline, controls, convergence ==")
    for name, (alpha, T) in CASES.items():
        d = GeneralDrive(S=Z, alpha_ac=AAC, amp=ramp(T))
        o = run(d, alpha, T)
        print(f"  {name:4s} alpha={alpha} T={T}: {fmt(o)}")
    for name in ("FAST", "STEP"):
        alpha, T = CASES[name]
        d = GeneralDrive(S=Z, alpha_ac=AAC, amp=ramp(T))
        a, b = run(d, alpha, T, law=False), run(d, alpha, T, nz=2 * nz_for(d), law=False)
        print(f"  convergence {name}: nz {nz_for(d)} vs x2: I_max {a['I_max']:.5f}/{b['I_max']:.5f}, "
              f"I_fin {a['I_fin']:.2e}/{b['I_fin']:.2e}")
    # step to S below the zero: n = 1 bright at the end (what a non-forbidden final looks like)
    d = GeneralDrive(S=2.5, alpha_ac=AAC, amp=ramp(0.5))
    print(f"  reference S_f=2.5 (J1={jv(1, 2.5):.3f}), FAST: {fmt(run(d, 1.0, 0.5))}")


def section_atten():
    print("== 1. attenuation A -> A e^{-gamma zeta}; amplitude loss across L = 1 - e^{-gamma} ==")
    for name in ("FAST", "MID", "STEP", "SLOW"):
        alpha, T = CASES[name]
        for loss in (0.0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.4, 0.6):
            d = GeneralDrive(S=Z, alpha_ac=AAC, amp=ramp(T), gamma=-np.log(1 - loss))
            print(f"  {name:4s} loss {loss:4.0%}: {fmt(run(d, alpha, T))}")


def section_disp():
    print("== 2. dispersion: envelope broadened by Gaussian of width sig1*sqrt(zeta) ==")
    for name in ("FAST", "MID", "STEP"):
        alpha, T = CASES[name]
        for sig1 in (0.0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8):
            amp = dispersed_ramp(T, sig1) if sig1 else ramp(T)
            d = GeneralDrive(S=Z, alpha_ac=AAC, amp=amp)
            print(f"  {name:4s} sig1 {sig1:4.2f}: {fmt(run(d, alpha, T, extra_t=0.5 + 3 * sig1))}")


def section_nl():
    print("== 3. nonlinearity: second harmonic m2 = eps2 * zeta * S * A^2, phase psi2 ==")
    for name in ("FAST", "MID", "STEP"):
        alpha, T = CASES[name]
        for eps2 in (0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2):
            d = GeneralDrive(S=Z, alpha_ac=AAC, amp=ramp(T), eps2=eps2)
            print(f"  {name:4s} eps2 {eps2:5.3f} psi2 0   : {fmt(run(d, alpha, T))}")
        for eps2 in (0.02, 0.1):
            d = GeneralDrive(S=Z, alpha_ac=AAC, amp=ramp(T), eps2=eps2, psi2=np.pi / 2)
            print(f"  {name:4s} eps2 {eps2:5.3f} psi2 pi/2: {fmt(run(d, alpha, T, law=False))}")


def section_ring():
    print("== 4. transducer ringing after the ramp: R e^{-x/tau} sin(2 pi f x + psi), x = s - T ==")
    for name in ("FAST", "MID", "STEP"):
        alpha, T = CASES[name]
        for R, tau, f, psi in [(0.0, 1, 1, 0), (0.02, 0.2, 3, 0), (0.05, 0.2, 3, 0), (0.1, 0.2, 3, 0),
                               (0.2, 0.2, 3, 0), (0.1, 0.05, 3, 0), (0.1, 1.0, 3, 0), (0.1, 3.0, 3, 0),
                               (0.1, 0.2, 0.5, 0), (0.1, 0.2, 10, 0), (0.1, 0.2, 3, np.pi / 2),
                               (0.1, 0.2, 3, np.pi)]:
            amp = ringing_ramp(T, R, tau, f, psi) if R else ramp(T)
            d = GeneralDrive(S=Z, alpha_ac=AAC, amp=amp)
            o_meas = run(d, alpha, T, extra_t=0.5)                      # fixed measurement window
            o_late = run(d, alpha, T, extra_t=0.5 + 6 * tau, law=False)  # after the ringing has died
            print(f"  {name:4s} R {R:.2f} tau {tau:4.2f} f {f:4.1f} psi {psi:4.2f}: {fmt(o_meas)} "
                  f"|| late I_fin {o_late['I_fin']:.2e}, C {o_late['C']:.1f}")


def section_combo():
    print("== 5. combined ==")
    for name in ("FAST", "MID"):
        alpha, T = CASES[name]
        for label, kw in [("atten 10% + ringing R .05 tau .2 f 3",
                           dict(amp=ringing_ramp(T, 0.05, 0.2, 3.0), gamma=-np.log(0.9))),
                          ("dispersion sig1 .1 + NL eps2 .02",
                           dict(amp=dispersed_ramp(T, 0.1), eps2=0.02)),
                          ("atten 5% + NL eps2 .02 + ringing R .02",
                           dict(amp=ringing_ramp(T, 0.02, 0.2, 3.0), gamma=-np.log(0.95), eps2=0.02))]:
            d = GeneralDrive(S=Z, alpha_ac=AAC, **kw)
            print(f"  {name:4s} {label}: {fmt(run(d, alpha, T))}")


def section_retune():
    print("== 6. essence test: with attenuation, retune S so that int_0^1 J1(S e^{-gamma z}) dz = 0 ==")
    z = np.linspace(0, 1, 4001)
    for loss in (0.1, 0.2, 0.4):
        g = -np.log(1 - loss)
        Sst = brentq(lambda S: np.trapezoid(jv(1, S * np.exp(-g * z)), z), 3.5, 6.0)
        frac = np.mean(np.abs(jv(1, Sst * np.exp(-g * z))) > 0.05)
        for name in ("FAST", "STEP"):
            alpha, T = CASES[name]
            d0 = GeneralDrive(S=Z, alpha_ac=AAC, amp=ramp(T), gamma=g)
            d1 = GeneralDrive(S=Sst, alpha_ac=AAC, amp=ramp(T), gamma=g)
            print(f"  loss {loss:.0%} {name}: S=j11 -> {fmt(run(d0, alpha, T))}")
            print(f"  loss {loss:.0%} {name}: S*={Sst:.4f} (|J1|>0.05 on {frac:.0%} of depth) -> "
                  f"{fmt(run(d1, alpha, T))}")


def section_integ():
    print("== 7. integrated satellite (window |eta - q| < q/2) vs centre ==")
    q = 2 * np.pi * AAC
    eta = np.linspace(0.5 * q, 1.5 * q, 801)
    alpha, T = CASES["FAST"]
    for label, kw in [("ideal", {}), ("atten 20%", dict(gamma=-np.log(0.8))),
                      ("NL eps2 .05", dict(eps2=0.05)), ("ring R .1 tau .2", dict(amp=ringing_ramp(T, 0.1, 0.2, 3.0)))]:
        kw.setdefault("amp", ramp(T))
        d = GeneralDrive(S=Z, alpha_ac=AAC, **kw)
        oc, oi = run(d, alpha, T), run(d, alpha, T, eta=eta, law=False)
        print(f"  FAST {label:16s}: centre C {oc['C']:9.1f} | integrated I_max {oi['I_max']:.3f} "
              f"I_fin {oi['I_fin']:.2e} C {oi['C']:9.1f}")


if __name__ == "__main__":
    secs = sys.argv[1:] or ["base", "atten", "disp", "nl", "ring", "combo", "retune", "integ"]
    for s in secs:
        globals()[f"section_{s}"]()
