"""Two-beam Takagi-Taupin for one frozen phase profile -> I(eta).

Everything is dimensionless:

- depth      ``zeta = z / L`` in [0, 1]
- coupling   ``alpha = L / Lambda_ext``  (Pendelloesung periods across the slab)
- rocking    ``eta``, the deviation parameter already multiplied by L

The common forward susceptibility ``chi_0`` is gauged out exactly: in symmetric
geometry it enters both equations with the same coefficient, so a shared phase
substitution removes it. Complex ``chi_0`` only leaves a common ``exp(-mu z)``,
which cancels in every normalized observable.

Absorption. In Laue a complex ``alpha`` carries the Borrmann (anomalous) part; in
Bragg absorption: normal (chi_0'') enters as ``eta -> eta + i kappa``, anomalous (chi_h'',
Borrmann) as ``c -> c + kappa_h / 2`` (see ``solve_bragg``).
A real ``alpha`` and ``kappa = 0`` are the non-absorbing case. See DESIGN.md for
the full assumption list.

Cost. Each RK4 step needs the phase factor ``e^{i phi}`` at the start, middle and
end of the step. The end of one step is the start of the next, so the factor is
evaluated twice per step, not four times.

Laue (transmission). Both beams travel +z; entrance ``D0(0)=1, Dh(0)=0``::

    dD0/dzeta = i pi alpha e^{+i phi} Dh
    dDh/dzeta = i pi alpha e^{-i phi} D0 + i eta Dh

Bragg (reflection). ``Dh`` travels -z, so ``d/ds_h = -d/dzeta``, and the boundary
conditions sit at opposite faces (``D0(0)=1``, ``Dh(1)=0``) -- a two-point BVP.
Rather than shoot, integrate the Riccati variable ``R = Dh / D0`` backward from
``R(1) = 0``, which turns it into an IVP in the stable direction::

    dR/dzeta = -i pi alpha e^{-i phi} - i eta R - i pi alpha e^{+i phi} R^2

Analytic ground truth, asserted in ``tests/test_core.py``:

- Laue, ``phi=0``, ``eta=0``  ->  ``|Dh(1)|^2 = sin^2(pi alpha)``
- Laue, real ``eta``          ->  ``|D0|^2 + |Dh|^2 = 1`` (lossless)
- Bragg, ``phi=0``            ->  total-reflection plateau for ``|eta| < 2 pi alpha``
"""

from __future__ import annotations

from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "PhaseField",
    "ZERO_PHASE",
    "required_nz",
    "solve_laue",
    "solve_bragg",
    "amplitude",
    "rocking",
]


class PhaseField(Protocol):
    """Frozen phase profile ``phi(zeta)``, vectorized over depth."""

    def __call__(self, zeta: NDArray[np.float64]) -> NDArray[np.float64]: ...


def ZERO_PHASE(zeta):
    """Unmodulated reference."""
    return np.zeros_like(np.asarray(zeta, dtype=np.float64))


def required_nz(
    alpha: float,
    alpha_ac: float,
    per_cycle: int = 20,
    S: float = 0.0,
    eta_max: float = 0.0,
) -> int:
    """Depth steps needed to resolve the fastest oscillation in the integrand.

    The grating turns over ``alpha_ac`` times across the slab and the
    Pendelloesung solution ``alpha`` times. But the integrand also carries
    ``e^{-i phi}``, whose phase sweeps up to ``S q`` per unit depth, and
    ``e^{i eta zeta}``, which turns ``eta / 2 pi`` times; at large ``S`` and
    high channels these dominate. Pass ``S`` and the largest ``|eta|`` to
    include them: e.g. ``S = 6`` read at ``eta = 5 q`` needs ~4x the steps the
    ``alpha_ac`` rule gives (notebook Entry 31).
    """
    cycles = max(alpha, alpha_ac, eta_max / (2.0 * np.pi) + S * alpha_ac)
    return max(256, int(per_cycle * cycles))


def _state_shape(eta_a: NDArray[np.float64], phase) -> tuple[int, ...]:
    """Broadcast shape of the integrator state.

    A phase profile that returns a scalar per depth gives the plain ``eta`` shape.
    A profile that returns an extra leading axis -- e.g. ``Drive.phi`` handed a
    column of times ``t[:, None]`` -- gives ``(n_t,) + eta.shape``, so one RK4
    sweep produces a whole time series instead of one snapshot. That matters: the
    per-step cost is dominated by Python loop overhead, not by array width.
    """
    probe = np.asarray(phase(np.float64(0.0)), dtype=np.float64)
    return np.broadcast_shapes(eta_a.shape, probe.shape)


def _couplings(phi, c, zeta: float) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """``(c e^{+i phi}, c e^{-i phi})`` at one depth."""
    e = np.exp(1j * np.asarray(phi(np.float64(zeta)), dtype=np.float64))
    return c * e, c * np.conj(e)


def _check(alpha, n_z: int) -> None:
    if np.real(alpha) < 0:
        raise ValueError(f"alpha must be non-negative, got {alpha}")
    if n_z < 1:
        raise ValueError(f"n_z must be >= 1, got {n_z}")


def solve_laue(
    alpha: complex,
    eta: ArrayLike = 0.0,
    phi: PhaseField | None = None,
    n_z: int = 512,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Integrate the Laue system from ``zeta=0`` to ``1`` (fixed-step RK4).

    Vectorized over ``eta``: the phase profile is shared, so a whole rocking
    curve costs one sweep. Returns the exit amplitudes ``(D0, Dh)``.

    A complex ``alpha`` models the Borrmann part of absorption in a
    centrosymmetric crystal (``chi_h = chi_hbar``, so the same coupling enters
    both equations). The common ``exp(-mu z)`` of ``chi_0''`` is still gauged out.
    """
    _check(alpha, n_z)
    eta_a = np.asarray(eta, dtype=np.float64)
    phase = ZERO_PHASE if phi is None else phi
    shape = _state_shape(eta_a, phase)
    D0 = np.ones(shape, dtype=np.complex128)
    Dh = np.zeros(shape, dtype=np.complex128)
    c = 1j * np.pi * alpha
    ie = 1j * eta_a
    h = 1.0 / n_z

    def f(cp, cm, d0, dh):
        return cp * dh, cm * d0 + ie * dh

    start = _couplings(phase, c, 0.0)
    for i in range(n_z):
        mid, end = _couplings(phase, c, (i + 0.5) * h), _couplings(phase, c, (i + 1) * h)
        k1a, k1b = f(*start, D0, Dh)
        k2a, k2b = f(*mid, D0 + 0.5 * h * k1a, Dh + 0.5 * h * k1b)
        k3a, k3b = f(*mid, D0 + 0.5 * h * k2a, Dh + 0.5 * h * k2b)
        k4a, k4b = f(*end, D0 + h * k3a, Dh + h * k3b)
        D0 = D0 + (h / 6.0) * (k1a + 2 * k2a + 2 * k3a + k4a)
        Dh = Dh + (h / 6.0) * (k1b + 2 * k2b + 2 * k3b + k4b)
        start = end

    return D0, Dh


def solve_bragg(
    alpha: float,
    eta: ArrayLike = 0.0,
    phi: PhaseField | None = None,
    n_z: int = 512,
    kappa: float = 0.0,
    kappa_h: float = 0.0,
) -> NDArray[np.complex128]:
    """Integrate the Bragg Riccati equation backward from ``R(1)=0`` to ``zeta=0``.

    Returns the complex amplitude reflectivity ``R(0) = Dh(0)/D0(0)`` at the
    entrance surface; the measured reflectivity is ``|R|^2``.

    ``kappa >= 0`` is normal absorption (chi_0''), which in Bragg does *not* gauge
    out: the two beams travel in opposite directions, so it enters ``R`` as
    ``eta -> eta + i kappa``. A layer at depth ``zeta`` then contributes with weight
    ``exp(-kappa zeta)``, i.e. ``kappa`` is the round-trip amplitude attenuation per
    unit depth (``mu L (1/gamma_0 + 1/gamma_h) / 2``).

    ``kappa_h >= 0`` is anomalous absorption (chi_h'', the Borrmann effect) in the
    same units, ``kappa_h / kappa = chi_ih / chi_i0``. For a centrosymmetric crystal
    ``chi_h = chi_hbar`` is one complex number, so both couplings become
    ``c = i pi alpha + kappa_h / 2``: the two Bloch modes then absorb at different
    rates, and the static reflectivity is higher on the ``eta < 0`` edge (the
    low-absorption edge; the low-angle side when ``eta`` grows with ``theta -
    theta_B``). Thick slab, kappa > 0 -> the Darwin-Prins root of
    ``c R^2 + i (eta + i kappa) R + c = 0`` with this ``c`` (tests/test_core.py).
    Laue has no ``kappa_h`` (notebook Entry 77; DESIGN.md item 4).
    """
    _check(alpha, n_z)
    if kappa < 0:
        raise ValueError(f"kappa must be non-negative, got {kappa}")
    if kappa_h < 0:
        raise ValueError(f"kappa_h must be non-negative, got {kappa_h}")

    eta_a = np.asarray(eta, dtype=np.float64)
    phase = ZERO_PHASE if phi is None else phi
    R = np.zeros(_state_shape(eta_a, phase), dtype=np.complex128)
    c = 1j * np.pi * alpha + 0.5 * kappa_h
    ie = 1j * (eta_a + 1j * kappa)
    h = -1.0 / n_z  # backward

    def f(cp, cm, r):
        return -cm - (ie + cp * r) * r

    start = _couplings(phase, c, 1.0)
    for i in range(n_z):
        mid, end = (
            _couplings(phase, c, 1.0 + (i + 0.5) * h),
            _couplings(phase, c, 1.0 + (i + 1) * h),
        )
        k1 = f(*start, R)
        k2 = f(*mid, R + 0.5 * h * k1)
        k3 = f(*mid, R + 0.5 * h * k2)
        k4 = f(*end, R + h * k3)
        R = R + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        start = end

    if not np.all(np.isfinite(R)):
        raise FloatingPointError("Bragg Riccati integration diverged; raise n_z or lower alpha")
    return R


def amplitude(
    alpha: float,
    eta: ArrayLike,
    phi: PhaseField | None = None,
    geometry: str = "laue",
    n_z: int = 512,
) -> NDArray[np.complex128]:
    """Complex diffracted amplitude whose squared modulus is the intensity.

    ``geometry='laue'`` gives ``Dh(L)``; ``'bragg'`` gives ``R(0)``.

    The phase matters for more than bookkeeping. Subtracting the unmodulated
    solution at the *amplitude* level, ``Dh(S) - Dh(0)``, cancels the part of the
    crystal the acoustic wave has not reached: in the kinematic limit the
    integrand becomes ``(e^{-i phi} - 1)``, which is identically zero ahead of the
    wavefront, so the difference is supported only on the filled depth. The same
    subtraction done on intensities does not cancel anything.
    """
    if geometry == "laue":
        return solve_laue(alpha, eta, phi=phi, n_z=n_z)[1]
    if geometry == "bragg":
        return solve_bragg(alpha, eta, phi=phi, n_z=n_z)
    raise ValueError(f"geometry must be 'laue' or 'bragg', got {geometry!r}")


def rocking(
    alpha: float,
    eta: ArrayLike,
    phi: PhaseField | None = None,
    geometry: str = "laue",
    n_z: int = 512,
) -> NDArray[np.float64]:
    """Diffracted intensity ``I(eta)`` for one frozen snapshot.

    ``geometry='laue'`` returns ``|Dh(L)|^2``; ``'bragg'`` returns ``|R(0)|^2``.
    """
    return np.abs(amplitude(alpha, eta, phi=phi, geometry=geometry, n_z=n_z)) ** 2
