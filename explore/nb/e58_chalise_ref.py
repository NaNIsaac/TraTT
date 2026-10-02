"""Entry 58: reference case for Chalise et al. (arXiv:2603.28071) — does the Fourier weight of R(t) at a
detuning eta reproduce the strain spectrum at the matching frequency?

Dimensionless symmetric Bragg (tratt.core.solve_bragg; non-absorbing, slab of thickness 1, alpha = L/Lambda).
A packet enters at the surface at t = 0 and travels at v = 1 (t in transit units):
    phi(zeta, t) = S g(t - zeta),   g(tau) = exp(-tau / tau_d) sin(q0 tau) for tau > 0, else 0,
so the displacement returns to zero behind the packet and the strain eps ~ g' starts with a step
(cf. their Fig. 2a). A component of wavenumber K = omega (v = 1) is Bragg-matched at eta = +-K, so
R(t) at detuning eta is expected to oscillate at |omega| = |eta|.

Weight W(eta) = |Fourier component of R(t; eta) at omega = |eta|| (magnitude, as in the paper).
Compared with |eps^(omega)|, |phi^(omega)| = |eps^|/omega and |eps^|/omega^2 (each normalised to its peak),
to identify which normalisation makes the shapes coincide.
Reference: S << 1 and q0 >> Darwin half-width 2 pi alpha. Then one variable at a time:
  q0 sweep toward the Darwin scale (R_D = 2 pi alpha / q0) at fixed S.
Usage: python e58_chalise_ref.py [ref] [abs] [darwin]
"""
import os
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from tratt.core import required_nz, solve_bragg  # noqa: E402


def packet(S, q0, ncyc=3.0):
    tau_d = ncyc * 2 * np.pi / q0

    def g(tau):
        return np.where(tau > 0, np.exp(-np.clip(tau, 0, None) / tau_d) * np.sin(q0 * tau), 0.0)
    return g, tau_d


def run(alpha, S, q0, eta, T=0.6, nt=1200, nz=None, kappa=0.0):
    g, tau_d = packet(S, q0)
    t = np.linspace(0.0, T, nt, endpoint=False)
    nz = nz or required_nz(alpha, q0 / (2 * np.pi), per_cycle=24, S=S, eta_max=np.abs(eta).max())

    def phi(zeta):
        return S * g(t[:, None] - zeta)
    R = np.abs(solve_bragg(alpha, eta, phi=phi, n_z=nz, kappa=kappa)) ** 2          # (nt, n_eta)
    # Fourier component exactly at omega = |eta| (a nearest-FFT-bin pick adds ~15 % leakage noise, Entry 59)
    W = np.abs(np.mean((R - R.mean(axis=0)) * np.exp(1j * np.abs(eta)[None, :] * t[:, None]), axis=0))
    # strain and displacement spectra of the packet (time axis at fixed depth = space axis, v = 1)
    tau = np.linspace(0, 40 * tau_d, 20001)
    gg = g(tau)
    eps = np.gradient(gg, tau)
    om = np.abs(eta)
    E = np.abs([np.trapezoid(eps * np.exp(1j * o * tau), tau) for o in om])
    P = np.abs([np.trapezoid(gg * np.exp(1j * o * tau), tau) for o in om])
    return dict(eta=eta, W=W, eps=E, phi=P, eps_w2=E / om**2, nz=nz, R0=R.mean(axis=0))


def shape_err(a, b):
    a, b = a / a.max(), b / b.max()
    return float(np.max(np.abs(a - b)))


def report(lab, r):
    print(f" {lab}: n_z = {r['nz']}; peak of W at eta = {r['eta'][np.argmax(r['W'])]:.0f}")
    r["W_r0"] = r["W"] / np.sqrt(r["R0"])                    # heterodyne reference divided out
    for k in ("eps", "phi", "eps_w2"):
        print(f"   max |W/maxW - {k}/max{k}| = {shape_err(r['W'], r[k]):.3f}")


def section_ref():
    alpha, q0 = 2.0, 2 * np.pi * 40
    eta = np.linspace(40, 3 * q0, 160)
    out = {}
    for S in (0.02, 0.05):
        for sgn in (+1, -1):
            r = run(alpha, S, q0, sgn * eta)
            report(f"S={S} eta sign {sgn:+d} (R_D = 2 pi alpha / q0 = {2*np.pi*alpha/q0:.3f})", r)
            out[f"S{S}_{sgn}"] = {k: np.asarray(v).tolist() for k, v in r.items()}
    r2 = run(alpha, 0.02, q0, eta, nt=2400, nz=2 * required_nz(alpha, 40, per_cycle=24, S=0.02, eta_max=eta.max()))
    report("convergence (nt x2, n_z x2), S=0.02, +eta", r2)
    print(f"   W change vs base: {np.max(np.abs(r2['W'] - np.array(out['S0.02_1']['W']))) / max(r2['W']):.2e}")
    json.dump(out, open(HERE / "e58_ref.json", "w"))


def section_abs():
    """Entry 59: same comparison with absorption, so the back face is invisible (Entry 58 artefact).
    alpha = 8 (Darwin half-width 2 pi alpha = 50), kappa = 6 (kappa / 2 pi alpha = 0.12, back face e^-6),
    Darwin ratio 2 pi alpha / q0 = 0.05."""
    alpha, kappa = 8.0, 6.0
    q0 = 2 * np.pi * alpha / 0.05
    eta = np.linspace(0.15 * q0, 3 * q0, 160)
    out = {}
    for S in (0.02, 0.05):
        r = run(alpha, S, q0, eta, kappa=kappa)
        report(f"abs: S={S}, kappa={kappa}, q0={q0:.0f}", r)
        for k in ("eps", "phi", "eps_w2"):
            print(f"   max |W_r0/max - {k}/max| = {shape_err(r['W_r0'], r[k]):.3f}")
        out[f"S{S}"] = {k: np.asarray(v).tolist() for k, v in r.items()}
    json.dump(out, open(HERE / "e58_abs.json", "w"))


def section_darwin():
    """Entry 59, one variable: the Darwin ratio R_D = 2 pi alpha / q0 (0.05 is section abs).
    Deficit D(eta) = W eta^2 / |eps^| relative to its large-eta value, vs eta / (2 pi alpha)."""
    alpha, kappa, S = 8.0, 6.0, 0.02
    w_d = 2 * np.pi * alpha
    out = {}
    for RD in (0.2, 0.5):
        q0 = w_d / RD
        eta = np.linspace(0.3 * w_d, 3 * q0 + 2 * w_d, 160)
        r = run(alpha, S, q0, eta, kappa=kappa)
        D = r["W"] * eta ** 2 / r["eps"]
        D = D / np.median(D[eta > 2.5 * w_d]) if np.any(eta > 2.5 * w_d) else D / D[-1]
        out[f"RD{RD}"] = {k: np.asarray(v).tolist() for k, v in r.items()} | {"D": D.tolist(), "q0": q0}
        sel = [np.argmin(abs(eta / w_d - x)) for x in (0.4, 0.7, 1.0, 1.5, 2, 3, 5, 8)]
        print(f" R_D={RD}: q0={q0:.0f}, n_z={r['nz']}; D(eta/2pi alpha) = " +
              "  ".join(f"{eta[i] / w_d:.1f}:{D[i]:.2f}" for i in sel), flush=True)
    json.dump(out, open(HERE / "e58_darwin.json", "w"))


if __name__ == "__main__":
    for s in (sys.argv[1:] or ["ref"]):
        globals()[f"section_{s}"]()
