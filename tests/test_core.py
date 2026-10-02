"""The analytic ground truth of DESIGN.md section 0, as asserts.

These exist so the solver never becomes a black box: every one of them is a
closed-form statement that must hold independently of the implementation.
"""

import numpy as np
import pytest
from scipy.special import jv

from tratt.core import required_nz, rocking, solve_bragg, solve_laue
from tratt.drive import Drive


# --- fact 1: Laue Pendelloesung -------------------------------------------------

@pytest.mark.parametrize("alpha", [0.25, 0.5, 0.9, 1.7])
def test_laue_pendelloesung(alpha):
    """Unmodulated, on Bragg: |Dh(1)|^2 = sin^2(pi alpha)."""
    I = rocking(alpha, 0.0, phi=None, geometry="laue", n_z=2048)
    assert I == pytest.approx(np.sin(np.pi * alpha) ** 2, abs=1e-8)


# --- fact 2: losslessness -------------------------------------------------------

def test_laue_is_lossless_with_modulation():
    """|D0|^2 + |Dh|^2 = 1 for real eta and any real phase profile."""
    drive = Drive(S=1.3, alpha_ac=5.0)
    eta = np.linspace(-40.0, 40.0, 33)
    D0, Dh = solve_laue(0.8, eta, phi=drive.phi(1.4), n_z=required_nz(0.8, 5.0))
    total = np.abs(D0) ** 2 + np.abs(Dh) ** 2
    assert np.allclose(total, 1.0, atol=1e-9)


# --- fact 3: satellites sit at analytically known eta ---------------------------

def test_satellites_land_on_predicted_channels():
    """Peaks of I(eta) coincide with drive.eta_channel(n) = -n*q."""
    alpha, S, alpha_ac, t = 0.05, 1.0, 4.0, 1.5
    drive = Drive(S=S, alpha_ac=alpha_ac)
    q = drive.q
    eta = np.linspace(-2.5 * q, 2.5 * q, 2001)
    I = rocking(alpha, eta, phi=drive.phi(t), geometry="laue", n_z=required_nz(alpha, alpha_ac))

    for n in (-2, -1, 0, 1, 2):
        centre = drive.eta_channel(n)
        window = np.abs(eta - centre) < 0.5 * q
        peak = eta[window][np.argmax(I[window])]
        # the local maximum must be within a tenth of the channel spacing
        assert abs(peak - centre) < 0.1 * q, f"channel {n}: peak {peak:.3f} vs {centre:.3f}"


# --- fact 4: kinematic limit reproduces the Bessel ladder -----------------------

def test_kinematic_limit_matches_the_born_integral():
    """For alpha -> 0 the solver must reproduce the single-scattering integral.

    ``Dh(1) = i pi alpha int_0^1 e^{-i phi(zeta)} e^{i eta (1-zeta)} dzeta``.
    This pins the whole solver+drive chain, including the sign of the eta term
    that fixes which side channel ``+n`` sits on.
    """
    alpha, S, alpha_ac, t = 1e-3, 0.9, 6.0, 1.5
    drive = Drive(S=S, alpha_ac=alpha_ac)
    eta = np.linspace(-3.0 * drive.q, 3.0 * drive.q, 1201)
    I = rocking(alpha, eta, phi=drive.phi(t), geometry="laue", n_z=required_nz(alpha, alpha_ac))

    zeta = np.linspace(0.0, 1.0, 20001)
    kernel = np.exp(-1j * drive.phi(t)(zeta))[None, :] * np.exp(1j * np.outer(eta, 1.0 - zeta))
    I_born = np.abs(1j * np.pi * alpha * np.trapezoid(kernel, zeta, axis=1)) ** 2

    assert np.max(np.abs(I - I_born)) < 1e-4 * I_born.max()


@pytest.mark.parametrize("S", [0.4, 0.9])
def test_kinematic_channel_heights_follow_the_bessel_ladder(S):
    """For alpha -> 0 the channel peak heights follow ``J_n(S)^2``.

    Read at the channel centers rather than integrated over a window: for integer
    ``alpha_ac`` the neighbouring sincs vanish there exactly, so this is leakage
    free.
    """
    alpha, alpha_ac, t = 1e-3, 6.0, 1.5
    drive = Drive(S=S, alpha_ac=alpha_ac)
    eta = np.array([drive.eta_channel(n) for n in (-2, -1, 0, 1, 2)])
    I = rocking(alpha, eta, phi=drive.phi(t), geometry="laue", n_z=required_nz(alpha, alpha_ac))
    I0 = I[2]

    for k, n in enumerate((-2, -1, 0, 1, 2)):
        want = (jv(n, S) / jv(0, S)) ** 2
        assert I[k] / I0 == pytest.approx(want, rel=0.01), f"n={n}: {I[k] / I0:.4g} vs {want:.4g}"


# --- fact 5: after transit the state is exactly periodic ------------------------

def test_step_onset_is_periodic_after_transit():
    """t >= 1 with a hard step: I(eta, t) repeats with the acoustic period.

    This is what makes t_f <= t_tr structural, so a sweep reporting t_f > t_tr
    under a step onset is a bug, not a discovery.
    """
    alpha, alpha_ac = 0.7, 3.0
    drive = Drive(S=1.1, alpha_ac=alpha_ac)
    eta = np.linspace(-2.0 * drive.q, 2.0 * drive.q, 201)
    n_z = required_nz(alpha, alpha_ac)

    t0 = 1.3
    a = rocking(alpha, eta, phi=drive.phi(t0), geometry="laue", n_z=n_z)
    b = rocking(alpha, eta, phi=drive.phi(t0 + drive.period), geometry="laue", n_z=n_z)
    assert np.allclose(a, b, atol=1e-10)

    # and it is genuinely still moving before the front exits
    c = rocking(alpha, eta, phi=drive.phi(0.4), geometry="laue", n_z=n_z)
    assert not np.allclose(a, c, atol=1e-3)


# --- Bragg: Darwin plateau ------------------------------------------------------

def test_bragg_total_reflection_plateau():
    """Unmodulated thick Bragg: |R|^2 -> 1 inside |eta| < 2 pi alpha, < 1 outside."""
    alpha = 4.0
    w = 2.0 * np.pi * alpha
    inside = np.linspace(-0.6 * w, 0.6 * w, 41)
    outside = np.array([-2.5 * w, -1.6 * w, 1.6 * w, 2.5 * w])

    R_in = np.abs(solve_bragg(alpha, inside, n_z=required_nz(alpha, alpha))) ** 2
    R_out = np.abs(solve_bragg(alpha, outside, n_z=required_nz(alpha, alpha))) ** 2

    assert np.all(R_in > 0.99)
    assert np.all(R_out < 0.5)
    assert np.all(R_in <= 1.0 + 1e-9)


def test_bragg_absorption_gives_darwin_prins_root():
    """Normal absorption kappa: an absorbing slab thick enough to hide its back face reflects like
    the semi-infinite crystal, the root |R| < 1 of c R^2 + i (eta + i kappa) R + c = 0."""
    alpha, kappa = 4.0, 12.0
    eta = np.linspace(-60.0, 60.0, 25)
    R = solve_bragg(alpha, eta, n_z=6000, kappa=kappa)
    c = 1j * np.pi * alpha
    roots = [np.roots([c, 1j * (e + 1j * kappa), c]) for e in eta]
    R_dp = np.array([r[np.argmin(np.abs(r))] for r in roots])
    assert np.max(np.abs(R - R_dp)) < 1e-5
    assert np.all(np.abs(R) < 1.0)          # absorption: no total reflection


def test_bragg_reflectivity_bounded_with_modulation():
    """Energy conservation: a lossless modulated Bragg crystal cannot reflect > 1."""
    drive = Drive(S=1.5, alpha_ac=4.0)
    eta = np.linspace(-3.0 * drive.q, 3.0 * drive.q, 301)
    R = np.abs(solve_bragg(2.0, eta, phi=drive.phi(1.7), n_z=required_nz(2.0, 4.0))) ** 2
    assert np.all(R <= 1.0 + 1e-9)


def test_bragg_fourier_weight_matches_first_order_transfer_function():
    """Entry 59: a weak strain packet entering an absorbing Bragg crystal makes R(t) oscillate at
    omega = |eta|; the Fourier weight equals the first-order closed form built on the static root R0
    (W = |conj(R0) A(|eta|) + R0 conj(A(-|eta|))| / T, A = i c S (1 - R0^2) G / (Lambda + i omega))."""
    alpha, kappa, S, q0 = 2.0, 5.0, 0.02, 2 * np.pi * 2.0 / 0.3
    c = 1j * np.pi * alpha
    tau_d = 3 * 2 * np.pi / q0

    def g(tau):
        return np.where(tau > 0, np.exp(-np.clip(tau, 0, None) / tau_d) * np.sin(q0 * tau), 0.0)

    T, nt = 2.4, 2400            # the window must outlast the slowest response near omega = q0
    t = np.linspace(0.0, T, nt, endpoint=False)
    eta = np.array([0.6, 1.0, 1.5, 2.5]) * q0
    R = np.abs(solve_bragg(alpha, eta, phi=lambda z: S * g(t[:, None] - z), n_z=1500, kappa=kappa)) ** 2
    W = np.abs(np.mean((R - R.mean(axis=0)) * np.exp(1j * eta[None, :] * t[:, None]), axis=0))

    tau = np.linspace(0, 40 * tau_d, 40001)

    def G(w):
        return np.trapezoid(g(tau) * np.exp(1j * w * tau), tau)

    pred = []
    for e in eta:
        r = np.roots([c, 1j * (e + 1j * kappa), c])
        R0 = r[np.argmin(np.abs(r))]
        lam = 1j * (e + 1j * kappa) + 2 * c * R0

        def A(w):
            return 1j * c * S * (1 - R0 ** 2) * G(w) / (lam + 1j * w)
        pred.append(abs(np.conj(R0) * A(e) + R0 * np.conj(A(-e))) / T)
    assert np.allclose(W / np.array(pred), 1.0, atol=0.03)
