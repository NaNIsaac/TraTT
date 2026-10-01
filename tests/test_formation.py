"""The steady reference and the formation-time machinery.

The load-bearing statement here is structural, not empirical: with a hard step
onset the traveling profile and the same-phase steady profile become *identical*
once the front has crossed the slab, so any sweep reporting t_f > t_tr at
rise = 0 is a bug in the pipeline, not a physical memory effect.
"""

import numpy as np
import pytest

from tratt.core import required_nz, rocking
from tratt.cube import front_and_full, rocking_cube
from tratt.drive import Drive
from tratt.observe import channel_stats, channel_window, formation_time


# --- the steady reference itself ------------------------------------------------

def test_steady_phi_is_periodic_in_time():
    d = Drive(S=0.8, alpha_ac=5.0, rise=0.3)
    z = np.linspace(0.0, 1.0, 501)
    a = d.steady_phi(0.37)(z)
    b = d.steady_phi(0.37 + d.period)(z)
    assert np.allclose(a, b, atol=1e-12)


def test_steady_phi_ignores_rise():
    z = np.linspace(0.0, 1.0, 201)
    a = Drive(S=0.8, alpha_ac=5.0, rise=0.0).steady_phi(0.4)(z)
    b = Drive(S=0.8, alpha_ac=5.0, rise=0.7).steady_phi(0.4)(z)
    assert np.allclose(a, b, atol=1e-15)


def test_front_vanishes_ahead_of_the_wavefront():
    d = Drive(S=0.5, alpha_ac=5.0)
    z = np.linspace(0.0, 1.0, 100001)
    t = 0.4
    assert np.max(np.abs(d.phi(t)(z[z > t]))) == 0.0
    assert np.max(np.abs(d.phi(t)(z[z < t]))) == pytest.approx(0.5, abs=1e-4)


@pytest.mark.parametrize("rise, t_bound", [(0.0, 1.0), (0.4, 1.4)])
def test_front_equals_steady_after_the_bound(rise, t_bound):
    """The structural upper bound on the settling time, at the profile level."""
    d = Drive(S=0.7, alpha_ac=4.0, rise=rise)
    z = np.linspace(0.0, 1.0, 2001)
    for t in (t_bound, t_bound + 0.13, t_bound + 1.0):
        assert np.max(np.abs(d.phi(t)(z) - d.steady_phi(t)(z))) < 1e-12
    # and genuinely different before it
    assert np.max(np.abs(d.phi(0.5 * t_bound)(z) - d.steady_phi(0.5 * t_bound)(z))) > 1e-3


def test_cubes_agree_after_transit():
    """Same bound, now through the solver: I_front = I_full for t >= 1."""
    alpha, d = 1.3, Drive(S=0.5, alpha_ac=5.0)
    eta = np.linspace(-1.5 * d.q, 1.5 * d.q, 201)
    times = np.array([1.0, 1.07, 1.5])
    If, Is = front_and_full(alpha, d, eta, times, n_z=required_nz(alpha, d.alpha_ac))
    assert np.max(np.abs(If - Is)) < 1e-12


# --- channel observables ---------------------------------------------------------

def test_channel_windows_tile_without_overlap():
    d = Drive(S=1.0, alpha_ac=5.0)
    eta = np.linspace(-2.5 * d.q, 2.5 * d.q, 1001)
    masks = [channel_window(eta, d, n) for n in (-2, -1, 0, 1, 2)]
    stacked = np.sum(masks, axis=0)
    assert stacked.max() == 1, "windows overlap"


def test_channel_stats_recover_a_known_gaussian():
    d = Drive(S=1.0, alpha_ac=5.0)
    eta = np.linspace(-0.5 * d.q, 0.5 * d.q, 4001)[1:-1]
    mu, sd = 3.0, 2.0
    I = (np.exp(-0.5 * ((eta - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi)))[None, :]
    c = channel_stats(eta, I, d, 0)
    assert c.p[0] == pytest.approx(1.0, rel=1e-6)
    assert c.eta_c[0] == pytest.approx(mu, abs=1e-6)
    assert c.sigma[0] == pytest.approx(sd, rel=1e-6)


def test_channel_stats_flag_an_empty_window():
    d = Drive(S=1.0, alpha_ac=5.0)
    eta = np.linspace(-0.5 * d.q, 0.5 * d.q, 101)
    c = channel_stats(eta, np.zeros((1, eta.size)), d, 0)
    assert c.p[0] == 0.0
    assert np.isnan(c.eta_c[0]) and np.isnan(c.sigma[0])


# --- formation time verdicts -----------------------------------------------------

def test_formation_time_requires_the_tolerance_to_hold_forever_after():
    t = np.linspace(0.0, 3.0, 601)
    p_full = np.ones_like(t)
    p_front = 1.0 + 0.4 * np.exp(-2.0 * t) * np.cos(9.0 * t)
    F = formation_time(t, p_front, p_full, period=0.5, eps=0.05)
    assert F.status == "formed"
    assert np.all(F.e[t >= F.t_f] <= 0.05 + 1e-12)
    assert F.e[t < F.t_f].max() > 0.05


def test_formation_time_reports_non_convergence_instead_of_the_last_sample():
    t = np.linspace(0.0, 3.0, 61)
    p_full = np.ones_like(t)
    F = formation_time(t, p_full + 0.5, p_full, period=0.5)
    assert F.status == "not-converged" and np.isnan(F.t_f)


def test_step_onset_formation_time_never_exceeds_the_transit():
    """End to end: the pipeline must respect the structural bound."""
    alpha, d = 1.3, Drive(S=0.5, alpha_ac=5.0)
    eta = np.linspace(-2.0 * d.q, 2.0 * d.q, 801)
    times = d.times(t_max=1.6, per_period=16)
    If, Is = front_and_full(alpha, d, eta, times, n_z=required_nz(alpha, d.alpha_ac))
    for n in (0, 1, -1):
        F = formation_time(times, channel_stats(eta, If, d, n).p,
                           channel_stats(eta, Is, d, n).p, d.period, eps=0.05)
        assert F.status == "formed", (n, F.status)
        assert F.t_f <= 1.0 + 1e-12, f"channel {n}: t_f = {F.t_f} exceeds the transit time"


# --- satellite formation = Laue thickness curve of the filled depth -------------

def test_satellite_formation_follows_filled_depth_pendelloesung():
    """Channel n inside the filled depth is a two-beam Laue slab of coupling
    alpha |J_n(S)| and thickness t, so the raw window integral must follow
    pi^2 alpha_n int_0^{2 pi alpha_n t} J_0 with no fit parameter, and the centre
    must follow sin^2(pi alpha_n t). Resolved channels only (2 pi alpha << q)."""
    from scipy.integrate import quad
    from scipy.special import jv

    S, alpha, n = 1.5, 2.0, 1
    d = Drive(S=S, alpha_ac=20.0)
    q = d.q
    eta = np.linspace((n - 0.5) * q, (n + 0.5) * q, 601)
    t = np.array([0.2, 0.45, 0.75])
    I = rocking_cube(alpha, d, eta, t, n_z=800)
    an = alpha * abs(jv(n, S))
    p = np.trapezoid(I, eta, axis=1)
    p_th = [np.pi**2 * an * quad(lambda x: jv(0, x), 0, 2 * np.pi * an * x)[0] for x in t]
    assert np.allclose(p, p_th, rtol=0.05), (p, p_th)
    centre = I[:, np.argmin(np.abs(eta - n * q))]
    assert np.allclose(centre, np.sin(np.pi * an * t) ** 2, atol=0.05), centre


def test_accumulated_coupling_law_all_channels_with_rise():
    """Centre intensity is sin^2 of the coupling accumulated along the depth,
    Theta_n = |int_0^1 J_n(S A(t - zeta)) dzeta|. Covers the n = 0 two-layer law
    (unfilled layer contributes J_0(0) = 1) and a finite rise time."""
    from scipy.special import jv

    S, alpha = 1.5, 2.0
    z = np.linspace(0.0, 1.0, 2001)
    for rise in (0.0, 0.4):
        d = Drive(S=S, alpha_ac=20.0, rise=rise)
        t = np.array([0.2037, 0.5037, 0.8037, 1.1037])
        I = rocking_cube(alpha, d, np.array([0.0, d.q, 2 * d.q]), t, n_z=1200)
        for j, n in enumerate((0, 1, 2)):
            Th = [abs(np.trapezoid(jv(n, S * d.envelope(x - z)), z)) for x in t]
            assert np.allclose(I[:, j], np.sin(np.pi * alpha * np.array(Th)) ** 2, atol=0.06), (rise, n)


def test_bragg_satellite_saturates_as_filled_depth_slab():
    """Bragg: the filled layer is a Bragg slab of coupling pi alpha |J_n(S)|."""
    from scipy.special import jv

    S, alpha = 1.5, 1.0
    d = Drive(S=S, alpha_ac=20.0)
    t = np.array([0.1037, 0.4037, 0.7037])
    I = rocking_cube(alpha, d, np.array([d.q]), t, geometry="bragg", n_z=1200)[:, 0]
    assert np.allclose(I, np.tanh(np.pi * alpha * abs(jv(1, S)) * t) ** 2, atol=0.02), I


def test_dark_centre_keeps_integrated_satellite():
    """At alpha_n t = 1 the centre is dark but the window integral is still
    pi^2 alpha_n int_0^{2 pi} J_0 (about half its maximum): the intensity is
    redistributed into side lobes, not lost."""
    from scipy.integrate import quad
    from scipy.special import jv

    S, alpha = 1.5, 2.0
    d = Drive(S=S, alpha_ac=40.0)
    an = alpha * abs(jv(1, S))
    eta = np.linspace(0.5 * d.q, 1.5 * d.q, 2001)
    I = rocking_cube(alpha, d, eta, np.array([1.0 / an]), n_z=2400)[0]
    centre = I[np.argmin(np.abs(eta - d.q))]
    p = np.trapezoid(I, eta)
    p_th = np.pi**2 * an * quad(lambda x: jv(0, x), 0, 2 * np.pi)[0]
    assert centre < 0.01, centre
    assert abs(p / p_th - 1) < 0.05, (p, p_th)
