"""One parameter point -> the time-resolved rocking curve ``I(eta, t)``.

Two cubes come out of every condition and they must share the eta grid, the
depth grid and the geometry, otherwise the residual between them mixes physics
with discretization error:

``front``  the real traveling front, ``phi = S A(t - zeta) sin(q zeta - q t)``
``full``   the same-phase steady reference, ``phi = S sin(q zeta - q t)``

Both are stroboscopic: the deformation is frozen at each t and the static
two-beam problem is solved in depth. There is no time derivative of the X-ray
field and no memory of earlier snapshots, so "formation" here means the approach
of the frozen-snapshot response to its periodic form, nothing more.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from tratt.core import amplitude, required_nz, rocking
from tratt.drive import Drive

__all__ = ["rocking_cube", "front_and_full", "satellite_cube"]


def _default_nz(alpha: float, drive: Drive, eta: NDArray[np.float64], n_z: int | None) -> int:
    return n_z or required_nz(alpha, drive.alpha_ac, S=drive.S, eta_max=float(np.max(np.abs(eta))))


def _over_time(solve, drive: Drive, eta, times, steady: bool) -> NDArray:
    """Stack ``solve(phi)`` over time, a block of times per RK4 sweep -> ``(n_times, n_eta)``.

    Handing the drive a column of times makes phi(zeta) return shape (n_t, 1),
    which broadcasts against eta, so one sweep covers several times at once.
    Doing the whole series in one sweep is *slower*, not faster: the RK4 inner
    loop is memory bound, and its temporaries fall out of cache. Measured on a
    1501-point eta grid, the cost per time sample bottoms out near 5e4 live
    elements, so the block is sized to that instead of to n_t.
    """
    # ponytail: the 5e4 figure is one machine's cache, not a law. If a sweep runs
    # on other hardware and speed matters, re-measure by timing a few block sizes.
    phase = drive.steady_phi if steady else drive.phi
    block = max(1, 50_000 // max(eta.size, 1))
    parts = []
    for i in range(0, times.size, block):
        tb = times[i : i + block]
        parts.append(np.reshape(solve(phase(tb[:, None])), (tb.size, eta.size)))
    return np.concatenate(parts, axis=0)


def rocking_cube(
    alpha: float,
    drive: Drive,
    eta: NDArray[np.float64],
    times: NDArray[np.float64],
    geometry: str = "laue",
    n_z: int | None = None,
    steady: bool = False,
) -> NDArray[np.float64]:
    """``I(eta, t)`` as an ``(n_times, n_eta)`` array.

    ``steady=True`` swaps the traveling front for the same-phase full-profile
    reference. Everything else is held fixed.
    """
    eta = np.asarray(eta, dtype=np.float64)
    n_z = _default_nz(alpha, drive, eta, n_z)
    return _over_time(lambda phi: rocking(alpha, eta, phi=phi, geometry=geometry, n_z=n_z),
                      drive, eta, np.asarray(times, dtype=np.float64), steady)


def front_and_full(
    alpha: float,
    drive: Drive,
    eta: NDArray[np.float64],
    times: NDArray[np.float64],
    geometry: str = "laue",
    n_z: int | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The two cubes a formation-time analysis needs, on identical grids."""
    n_z = _default_nz(alpha, drive, np.asarray(eta, dtype=np.float64), n_z)
    return (
        rocking_cube(alpha, drive, eta, times, geometry=geometry, n_z=n_z, steady=False),
        rocking_cube(alpha, drive, eta, times, geometry=geometry, n_z=n_z, steady=True),
    )


def satellite_cube(
    alpha: float,
    drive: Drive,
    eta: NDArray[np.float64],
    times: NDArray[np.float64],
    geometry: str = "laue",
    n_z: int | None = None,
    steady: bool = False,
) -> NDArray[np.float64]:
    """``|D_h(S) - D_h(0)|^2``: the diffraction the acoustic wave is responsible for.

    A channel window is not a satellite. The crystal keeps its full thickness at
    every instant, so the part the wavefront has not reached goes on producing
    full-thickness fringes, and those land inside the satellite windows: at
    alpha_ac = 8 the n = +-1 window is mostly them at early times. Measuring a
    width there returns the window width or a fringe width, not the satellite.

    Subtracting the unmodulated amplitude removes exactly that contribution. In
    the kinematic limit the integrand of the difference carries a factor
    ``e^{-i phi} - 1``, which vanishes identically ahead of the wavefront, so the
    difference is supported on the filled depth ``zeta < t`` alone and the
    satellite it leaves has the finite-size width of that depth.

    This is a constructed observable, not the measured intensity: an experiment
    sees ``|D_h(S)|^2``. Report it as the acoustically induced amplitude
    difference, and read widths off it rather than off the raw window.
    """
    eta = np.asarray(eta, dtype=np.float64)
    n_z = _default_nz(alpha, drive, eta, n_z)
    A0 = amplitude(alpha, eta, phi=None, geometry=geometry, n_z=n_z)
    return _over_time(lambda phi: np.abs(amplitude(alpha, eta, phi=phi, geometry=geometry, n_z=n_z) - A0) ** 2,
                      drive, eta, np.asarray(times, dtype=np.float64), steady)
