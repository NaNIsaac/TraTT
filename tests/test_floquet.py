"""The exact harmonic (Floquet) solver against full TT, and its internal consistency."""

import numpy as np
import pytest

from tratt.core import required_nz, solve_laue
from tratt.drive import Drive
from tratt.floquet import FloquetLaue, auto_K


@pytest.mark.parametrize("alpha,aac,S,t", [(1.0, 10.0, 1.5, 0.0), (2.0, 2.0, 3.0, 0.3), (3.0, 0.75, 2.0, 0.1)])
def test_floquet_matches_full_tt(alpha, aac, S, t):
    q = 2 * np.pi * aac
    K = auto_K(S) + int(2 * alpha / aac) + 2
    eta = np.linspace(-(S + 2) * q - 3 * alpha, (S + 2) * q + 3 * alpha, 801)
    fl = FloquetLaue(alpha, aac, S, K, t=t)
    D0, Dh = solve_laue(alpha, eta, phi=Drive(S=S, alpha_ac=aac).steady_phi(t),
                        n_z=2 * required_nz(alpha, aac, S=S, eta_max=np.abs(eta).max()))
    D0f, Dhf = fl.fields(eta)
    assert np.allclose(np.abs(Dhf) ** 2, np.abs(Dh) ** 2, atol=5e-5)
    assert np.allclose(np.abs(D0f) ** 2 + np.abs(Dhf) ** 2, 1.0, atol=1e-6)   # flux conservation


def test_floquet_truncation_converges_and_mask_removes_coupling():
    fl = FloquetLaue(2.0, 3.0, 2.5, auto_K(2.5) + 2)
    fl2 = FloquetLaue(2.0, 3.0, 2.5, auto_K(2.5) + 6)
    eta = np.linspace(-60, 60, 201)
    assert np.allclose(fl.rocking(eta), fl2.rocking(eta), atol=1e-8)
    off = FloquetLaue(2.0, 3.0, 2.5, 6, mask=lambda k, m: False)
    assert np.allclose(off.rocking(eta), 0.0, atol=1e-14)


def test_floquet_reduces_to_unmodulated_laue_at_S0():
    fl = FloquetLaue(1.3, 4.0, 0.0, 3)
    eta = np.linspace(-30, 30, 301)
    s = np.sqrt((np.pi * 1.3) ** 2 + eta**2 / 4)
    exact = (np.pi * 1.3) ** 2 * np.sin(s) ** 2 / s**2
    assert np.allclose(fl.rocking(eta), exact, atol=1e-10)


def test_bessel_zero_resurrection_third_order():
    """At S = j_{1,1} the n = 1 channel population follows sin^2|g_eff| with
    g_eff = pi alpha (R/2)^2 C_1(S), the lowest indirect path a0 -> b_m -> a_k -> b_-1
    (notebook Entry 45). Checked in the perturbative regime R = 0.1."""
    import sys
    from pathlib import Path

    from scipy.special import jn_zeros

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "explore" / "nb"))
    from e45_zero import C_n, channel_population

    S = jn_zeros(1, 1)[0]
    alpha, aac = 1.0, 10.0
    P = channel_population(alpha, aac, S, 1)[0]
    g = np.pi * alpha * (alpha / (2 * aac)) ** 2 * C_n(1, S)
    assert abs(P / np.sin(g) ** 2 - 1) < 0.02, (P, np.sin(g) ** 2)


def test_two_band_gap_equals_independent_coupling_when_resolved():
    """Exact two-band Floquet structure (monodromy): unitary, and the gap at eta = n q is
    2 pi alpha |J_n(S)| in the resolved regime (notebook Entry 46)."""
    import sys
    from pathlib import Path

    from scipy.special import jv

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "explore" / "nb"))
    from e46_bands import bands

    aac, alpha, S = 5.0, 0.25, 1.5
    q = 2 * np.pi * aac
    eta = np.array([0.0, q, 2 * q])
    lam, gap, mix, det = bands(alpha, aac, S, eta)
    assert np.allclose(det, 1.0, atol=1e-5)
    assert np.allclose(gap, [2 * np.pi * alpha * abs(jv(n, S)) for n in range(3)], rtol=0.03)


def test_effective_2x2_hamiltonian_third_order():
    """Schrieffer-Wolff reduction of the exact Floquet H onto {a_0, b_-n} (notebook Entry 52):
    off-diagonal g = pi alpha [J_{-n}(1 - W_n (R/2)^2 / 2) - (R/2)^2 C_n] with error O(R^4), and the
    resonance shift delta* = -pi alpha R B_n. Measured as half the minimum splitting of the exact
    Floquet branches and its eta position."""
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "explore" / "nb"))
    from e52_heff import exact, predictions

    errs = []
    for R in (0.05, 0.1):
        alpha, aac, S, n = R * 10.0, 10.0, 1.5, 1
        ge, de = exact(alpha, aac, S, n)
        g1, g3C, g3, ds = predictions(alpha, aac, S, n)
        assert abs(ge - abs(g3)) < 1e-3 * abs(ge - abs(g1))      # full third order >> first order
        assert abs(de - ds) < 0.01 * abs(ds)
        errs.append(abs(ge - abs(g3)) / (np.pi * alpha))
    assert 3.5 < np.log(errs[1] / errs[0]) / np.log(2) < 4.5         # relative error ~ R^4


def test_exact_acoustic_phase_average_is_sum_of_floquet_populations():
    """Notebook Entry 53: b_m(1; phi) = e^{i m phi} b_m(1; 0) (phi = q t), so the uniform average over the
    acoustic phase of the exit intensity is sum_m |b_m(1)|^2; a 4-point average is not converged at
    non-integer alpha_ac while 16 points are (checked against direct TT)."""
    from scipy.special import jn_zeros

    alpha, aac, S = 1.0, 20.25, jn_zeros(1, 1)[0]
    q = 2 * np.pi * aac
    eta = np.array([q])
    K = auto_K(S) + 6
    b0 = FloquetLaue(alpha, aac, S, K).propagate(eta)[0, 2 * K + 1:]
    ms = np.arange(-K, K + 1)
    phi = 0.7
    b1 = FloquetLaue(alpha, aac, S, K, t=phi / q).propagate(eta)[0, 2 * K + 1:]
    assert np.allclose(b1, np.exp(1j * ms * phi) * b0, atol=1e-12)
    exact = np.sum(np.abs(b0) ** 2)
    nz = 2 * required_nz(alpha, aac, S=S, eta_max=q)
    tt = lambda N: np.mean([abs(solve_laue(alpha, eta, phi=Drive(S=S, alpha_ac=aac).steady_phi(k / (N * aac)),  # noqa: E731
                                           n_z=nz)[1][0]) ** 2 for k in range(N)])
    assert abs(tt(16) / exact - 1) < 1e-3
    assert abs(tt(4) / exact - 1) > 0.1
