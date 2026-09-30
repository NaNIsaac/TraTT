"""``I(eta, t)`` -> channel observables -> formation times.

Channel ``n`` is the window ``W_n = {(n - 1/2) q < eta < (n + 1/2) q}`` around the
analytically known centre ``eta_n = n q``. Inside it,

    p_n(t)      = int_{W_n} I(eta, t) d eta            (integrated intensity)
    eta_c,n(t)  = int_{W_n} eta I d eta / p_n          (first moment)
    sigma_n(t)^2= int_{W_n} (eta - eta_c)^2 I d eta / p_n

``p_n`` is an integrated intensity, not a probability and not an occupancy: it is
never normalized by the sum over channels. ``eta_n`` is where the window is
anchored, not a claim about where the true maximum lies at strong coupling.

Formation time. With ``p_n^front`` the traveling-front result and ``p_n^full`` the
same-phase steady reference (see cube.front_and_full),

    r_n(t) = p_n^front(t) - p_n^full(t)
    a_n    = period average of p_n^full          (scale = window content)
           | <p_n^full> - p_n^static |           (scale = induced change)
    e_n(t) = |r_n(t)| / a_n

The second scale exists because a channel window is not a pure satellite: the
unmodulated rocking curve decays only as 1/eta^2, so its tail can dominate the
outer windows. It cancels in r_n but not in the normalizer.

and ``t_f`` is the earliest time such that ``e_n <= eps`` holds at that sample
*and at every later one*. Two cases are reported rather than faked:

- ``e_n`` never exceeds eps          -> t_f = times[0], "no formation detected"
- ``e_n`` still above eps at the end -> t_f = nan, "not converged in the window"

Returning the last sample in the second case would silently manufacture a
formation time equal to the observation length.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from tratt.drive import Drive

__all__ = ["channel_window", "channel_stats", "ChannelStats", "Formation",
           "formation_time", "peak_fwhm"]


@dataclass(frozen=True)
class ChannelStats:
    """Per-time observables of one channel. Arrays have the shape of ``times``."""

    n: int
    p: NDArray[np.float64]
    eta_c: NDArray[np.float64]
    sigma: NDArray[np.float64]


@dataclass(frozen=True)
class Formation:
    """Formation-time verdict for one channel, with the evidence attached."""

    t_f: float
    status: str  # "formed" | "no-formation-detected" | "not-converged" | "channel-too-weak"
    eps: float
    a: float
    e: NDArray[np.float64]

    @property
    def ok(self) -> bool:
        return self.status in ("formed", "no-formation-detected")


def channel_window(eta: NDArray[np.float64], drive: Drive, n: int) -> NDArray[np.bool_]:
    """Boolean mask of ``|eta - n q| < q / 2``."""
    return np.abs(np.asarray(eta, dtype=np.float64) - drive.eta_channel(n)) < 0.5 * drive.q


def channel_stats(
    eta: NDArray[np.float64],
    I: NDArray[np.float64],
    drive: Drive,
    n: int,
) -> ChannelStats:
    """Integrated intensity, centroid and width of channel ``n``.

    ``I`` is ``(n_times, n_eta)``. Where ``p_n`` underflows, the centroid and the
    width are NaN rather than a ratio of two small numbers.
    """
    eta = np.asarray(eta, dtype=np.float64)
    I = np.atleast_2d(np.asarray(I, dtype=np.float64))
    mask = channel_window(eta, drive, n)
    if mask.sum() < 3:
        raise ValueError(f"channel {n} covers {mask.sum()} eta samples; refine the eta grid")

    e, w = eta[mask], I[:, mask]
    p = np.trapezoid(w, e, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        eta_c = np.trapezoid(e * w, e, axis=1) / p
        var = np.trapezoid((e - eta_c[:, None]) ** 2 * w, e, axis=1) / p
        bad = ~np.isfinite(p) | (p <= np.finfo(float).tiny)
        eta_c = np.where(bad, np.nan, eta_c)
        var = np.where(bad, np.nan, var)
    return ChannelStats(n=n, p=p, eta_c=eta_c, sigma=np.sqrt(np.clip(var, 0.0, None)))


def _period_average(times: NDArray[np.float64], y: NDArray[np.float64], period: float) -> float:
    """Mean of ``y`` over the final whole period; whole span if it is shorter."""
    span = times[-1] - times[0]
    lo = times[-1] - period if span >= period else times[0]
    m = times >= lo
    if m.sum() < 2:
        return float(np.mean(y))
    return float(np.trapezoid(y[m], times[m]) / (times[m][-1] - times[m][0]))


def formation_time(
    times: NDArray[np.float64],
    p_front: NDArray[np.float64],
    p_full: NDArray[np.float64],
    period: float,
    eps: float = 0.05,
    a_min: float = 0.0,
    p_static: float | None = None,
) -> Formation:
    """Earliest time after which ``|p_front - p_full| / a <= eps`` holds for good.

    The residual ``p_front - p_full`` is clean: both cubes carry the same static
    rocking-curve tail inside the window, so it cancels. The normalizer is what
    needs care. The window content itself, ``<p_full>``, is not the satellite
    strength -- at alpha = 1.3, alpha_ac = 5 the n = +-2 window holds more
    unmodulated tail than satellite, which inflates ``a`` and biases ``t_f`` low.

    Pass ``p_static``, the window integral of the S = 0 control, to normalize by
    the acoustically *induced* change ``|<p_full> - p_static|`` instead. The
    absolute value is required: modulation depletes the n = 0 window rather than
    filling it, so the induced change is negative there.

    Neither choice separates satellite from tail rigorously -- in a coherent
    calculation the two interfere and do not add. ``a`` is a stated scale for the
    residual, not a decomposition, so every reported ``t_f`` must say which one
    it used.

    ``a_min`` rejects channels whose scale is too weak for the ratio to mean
    anything; the default of 0 accepts everything and leaves that judgement to
    the caller, which must then look at ``Formation.a``.
    """
    times = np.asarray(times, dtype=np.float64)
    a = _period_average(times, np.asarray(p_full, dtype=np.float64), period)
    if p_static is not None:
        a = abs(a - float(p_static))
    if not np.isfinite(a) or a <= a_min:
        return Formation(np.nan, "channel-too-weak", eps, a, np.full(times.shape, np.nan))

    e = np.abs(np.asarray(p_front, dtype=np.float64) - np.asarray(p_full, dtype=np.float64)) / a
    over = np.flatnonzero(e > eps)
    if over.size == 0:
        return Formation(float(times[0]), "no-formation-detected", eps, a, e)
    if over[-1] == times.size - 1:
        return Formation(np.nan, "not-converged", eps, a, e)
    return Formation(float(times[over[-1] + 1]), "formed", eps, a, e)


def demo() -> None:
    """Self-check: the three verdicts, and the t >= 1 + rise structural bound."""
    t = np.linspace(0.0, 3.0, 301)

    # decaying transient that recrosses eps once -> t_f is the LAST crossing
    p_full = np.ones_like(t)
    p_front = 1.0 + 0.4 * np.exp(-2.0 * t) * np.cos(9.0 * t)
    f = formation_time(t, p_front, p_full, period=0.5, eps=0.05)
    assert f.status == "formed", f
    assert np.all(np.abs(p_front - p_full)[t >= f.t_f] <= 0.05 + 1e-12)
    assert np.any(np.abs(p_front - p_full)[t < f.t_f] > 0.05)

    assert formation_time(t, p_full, p_full, 0.5).status == "no-formation-detected"
    assert formation_time(t, p_full + 0.5, p_full, 0.5).status == "not-converged"
    assert formation_time(t, p_front, p_full, 0.5, a_min=2.0).status == "channel-too-weak"

    # a_n really is the period average of an oscillating steady level
    osc = 1.0 + 0.2 * np.sin(2 * np.pi * t / 0.5)
    assert abs(formation_time(t, osc, osc, period=0.5).a - 1.0) < 1e-3

    # induced scale: a window that is 90% static tail must not shrink the residual
    tail = 9.0
    F_win = formation_time(t, p_front + tail, p_full + tail, 0.5)
    F_ind = formation_time(t, p_front + tail, p_full + tail, 0.5, p_static=tail)
    assert abs(F_ind.a - 1.0) < 1e-9
    assert abs(F_win.a - 10.0) < 1e-9
    assert F_ind.t_f > F_win.t_f, (F_win.t_f, F_ind.t_f)

    print("observe.demo: all checks passed")


if __name__ == "__main__":
    demo()


def peak_fwhm(
    eta: NDArray[np.float64],
    y: NDArray[np.float64],
    centre: float,
    half_width: float,
) -> tuple[float, float, float]:
    """Full width at half maximum of the dominant peak inside one window.

    Returns ``(fwhm, peak_height, peak_position)``; ``fwhm`` is NaN when the curve
    does not fall to half maximum inside the window, which is the honest answer
    for a feature wider than the window rather than a clipped number.

    Measured on the highest point in the window, walking outward to the
    half-maximum crossing and interpolating linearly. On a raw rocking curve this
    picks up whatever fringe happens to be tallest; use it on the amplitude-
    subtracted curve of ``cube.satellite_cube``, where the unmodulated crystal's
    fringes are gone.
    """
    eta = np.asarray(eta, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    m = np.abs(eta - centre) < half_width
    if m.sum() < 5:
        raise ValueError("window holds too few eta samples; refine the grid")
    e, v = eta[m], y[m]
    k = int(np.argmax(v))
    top = v[k]
    if not np.isfinite(top) or top <= 0.0:
        return np.nan, float(top), np.nan
    lvl = 0.5 * top

    def crossing(step: int) -> float:
        i = k
        while 0 <= i + step < v.size and v[i + step] > lvl:
            i += step
        j = i + step
        if not (0 <= j < v.size):
            return np.nan  # never falls to half max inside the window
        f = (v[i] - lvl) / (v[i] - v[j])
        return float(e[i] + f * (e[j] - e[i]))

    lo, hi = crossing(-1), crossing(+1)
    return hi - lo, float(top), float(e[k])
