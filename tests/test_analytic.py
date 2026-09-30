"""Closed-form solutions of the unmodulated problem, used as solver ground truth.

These go beyond tests/test_core.py in two ways:

1. They check the *whole* rocking curve against an exact expression, not a
   single point or a limiting regime.
2. They check the RK4 convergence *order*, which is what actually certifies the
   integrator. A conserved quantity holding at one n_z proves nothing about
   accuracy; a fourth-order error decay does.

Derivations (handoff sections 6.1 and 6.3):

Laue, phi = 0. Substituting D0 = e^{i eta zeta/2} f, Dh = e^{i eta zeta/2} g
symmetrizes the system into f' = -i(eta/2) f + i kappa g, g' = i kappa f +
i(eta/2) g. That matrix squares to -(kappa^2 + eta^2/4) I, so with
nu = sqrt(kappa^2 + eta^2/4),

    D0(1) = e^{i eta/2} [cos nu - i (eta / 2 nu) sin nu]
    Dh(1) = i e^{i eta/2} (kappa / nu) sin nu

Bragg, phi = 0, eta = 0. R = i tanh(kappa (1 - zeta)) satisfies both
R' = -i kappa (1 + R^2) and R(1) = 0, hence |R(0)|^2 = tanh^2(pi alpha).
"""

import numpy as np
import pytest

from tratt.core import required_nz, rocking, solve_bragg, solve_laue
from tratt.drive import Drive


def laue_exact(alpha: float, eta: np.ndarray):
    """Exact (D0, Dh) at the exit face for the unmodulated Laue case."""
    kappa = np.pi * alpha
    eta = np.asarray(eta, dtype=np.float64)
    nu = np.sqrt(kappa**2 + 0.25 * eta**2)
    # nu -> 0 only when kappa = eta = 0; sinc keeps the limit finite there.
    sinc_nu = np.sinc(nu / np.pi)  # = sin(nu)/nu, continuous at nu = 0
    common = np.exp(0.5j * eta)
    D0 = common * (np.cos(nu) - 0.5j * eta * sinc_nu)
    Dh = 1j * common * kappa * sinc_nu
    return D0, Dh


# --- full rocking curve, not just the centre ------------------------------------

@pytest.mark.parametrize("alpha", [0.3, 1.0, 2.7])
def test_laue_full_rocking_curve_matches_closed_form(alpha):
    # n_z is set by the RK4 truncation error, which scales like h^4 eta^5:
    # at |eta| = 60 and n_z = 4096 it is ~3e-9, so 8192 buys the 1e-9 tolerance.
    eta = np.linspace(-60.0, 60.0, 601)
    D0, Dh = solve_laue(alpha, eta, phi=None, n_z=8192)
    D0_ex, Dh_ex = laue_exact(alpha, eta)
    assert np.max(np.abs(D0 - D0_ex)) < 1e-9
    assert np.max(np.abs(Dh - Dh_ex)) < 1e-9


@pytest.mark.parametrize("alpha", [0.25, 1.0, 4.0])
def test_bragg_centre_matches_tanh_squared(alpha):
    R = solve_bragg(alpha, 0.0, phi=None, n_z=4096)
    assert np.abs(R) ** 2 == pytest.approx(np.tanh(np.pi * alpha) ** 2, abs=1e-10)


# --- the integrator is genuinely fourth order -----------------------------------

def _laue_error(alpha, eta, n_z):
    D0, Dh = solve_laue(alpha, eta, phi=None, n_z=n_z)
    D0_ex, Dh_ex = laue_exact(alpha, eta)
    return max(np.max(np.abs(D0 - D0_ex)), np.max(np.abs(Dh - Dh_ex)))


def test_laue_error_decays_at_fourth_order():
    """Halving the step must cut the error by ~16. Anything less is a bug.

    Starts at n_z = 64: the 32 -> 64 ratio is still pre-asymptotic (~3.5).
    """
    alpha, eta = 1.7, np.linspace(-40.0, 40.0, 81)
    errs = [_laue_error(alpha, eta, n) for n in (64, 128, 256, 512)]
    orders = [np.log2(a / b) for a, b in zip(errs[:-1], errs[1:])]
    assert all(o > 3.8 for o in orders), f"errors {errs}, observed orders {orders}"


def test_losslessness_error_shrinks_with_the_grid():
    """Unitarity defect must *converge*, not merely be small at one n_z."""
    drive = Drive(S=1.3, alpha_ac=5.0)
    eta = np.linspace(-40.0, 40.0, 33)
    defects = []
    for n_z in (64, 128, 256):
        D0, Dh = solve_laue(0.8, eta, phi=drive.phi(1.4), n_z=n_z)
        defects.append(np.max(np.abs(np.abs(D0) ** 2 + np.abs(Dh) ** 2 - 1.0)))
    assert defects[0] > defects[1] > defects[2], defects


# --- structural invariants -------------------------------------------------------

def test_zero_coupling_is_transparent():
    eta = np.linspace(-20.0, 20.0, 51)
    D0, Dh = solve_laue(0.0, eta, phi=Drive(S=1.0, alpha_ac=3.0).phi(0.7), n_z=256)
    assert np.allclose(D0, 1.0, atol=1e-12)
    assert np.allclose(Dh, 0.0, atol=1e-12)


def test_constant_phase_offset_leaves_intensity_unchanged():
    """phi = phi0 rotates Dh only; with Dh(0) = 0 the exit intensity is invariant."""
    alpha, eta = 1.1, np.linspace(-30.0, 30.0, 121)
    I0 = rocking(alpha, eta, phi=None, geometry="laue", n_z=1024)
    I1 = rocking(alpha, eta, phi=lambda z: np.full_like(np.asarray(z, float), 0.9),
                 geometry="laue", n_z=1024)
    assert np.allclose(I0, I1, atol=1e-10)
