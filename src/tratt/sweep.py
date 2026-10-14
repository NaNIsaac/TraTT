"""Parameter grid -> formation-time table.

One row per (alpha, S, channel). The raw cubes are not kept: a grid of them is
tens of GB, and every observable here is reproducible from the recorded metadata
by rerunning ``condition``.

Each condition needs three solves on one shared eta grid and one shared depth
grid: the traveling front, the same-phase steady reference, and the S = 0 static
control that sets the induced normalizer. Mixing grids between them would put
discretization error into the residual.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import NDArray

from tratt.core import required_nz, rocking
from tratt.cube import front_and_full
from tratt.drive import Drive
from tratt.observe import channel_stats, formation_time

__all__ = ["Grid", "condition", "sweep"]


@dataclass(frozen=True)
class Grid:
    """Numerical resolution, recorded with every result so a row is reproducible."""

    n_eta: int = 1201
    eta_span: float = 3.0  # in units of q
    t_max: float = 2.0
    per_period: int = 32
    n_z: int | None = None  # None -> required_nz(alpha, alpha_ac, S=S, eta_max=eta_span*q)


def condition(
    alpha: float,
    S: float,
    alpha_ac: float = 5.0,
    rise: float = 0.0,
    geometry: str = "laue",
    channels: tuple[int, ...] = (0, 1, -1, 2, -2),
    eps: tuple[float, ...] = (0.02, 0.05, 0.10),
    grid: Grid = Grid(),
) -> list[dict]:
    """Formation times for one parameter point, one row per channel."""
    d = Drive(S=S, alpha_ac=alpha_ac, rise=rise)
    eta = np.linspace(-grid.eta_span * d.q, grid.eta_span * d.q, grid.n_eta)
    n_z = grid.n_z or required_nz(alpha, alpha_ac, S=S, eta_max=grid.eta_span * d.q)
    times = d.times(t_max=grid.t_max, per_period=grid.per_period)

    I_front, I_full = front_and_full(alpha, d, eta, times, geometry=geometry, n_z=n_z)
    I_static = rocking(alpha, eta, phi=None, geometry=geometry, n_z=n_z)[None, :]

    rows = []
    for n in channels:
        cf = channel_stats(eta, I_front, d, n)
        cs = channel_stats(eta, I_full, d, n)
        p_static = float(channel_stats(eta, I_static, d, n).p[0])
        row = dict(
            alpha=alpha, S=S, alpha_ac=alpha_ac, rise=rise, geometry=geometry, n=n,
            n_z=n_z, **{f"grid_{k}": v for k, v in asdict(grid).items() if k != "n_z"},
            p_static=p_static,
            p_steady_mean=float(np.mean(cs.p)),
            sigma_steady_mean=float(np.nanmean(cs.sigma)),
            eta_c_steady_mean=float(np.nanmean(cs.eta_c)),
        )
        for e in eps:
            F_ind = formation_time(times, cf.p, cs.p, d.period, eps=e, p_static=p_static)
            F_win = formation_time(times, cf.p, cs.p, d.period, eps=e)
            row[f"t_f_induced_{e:g}"] = F_ind.t_f
            row[f"status_induced_{e:g}"] = F_ind.status
            row[f"t_f_window_{e:g}"] = F_win.t_f
            row["a_induced"], row["a_window"] = F_ind.a, F_win.a   # the scales do not depend on eps
        rows.append(row)
    return rows


def sweep(
    alphas: NDArray[np.float64],
    Ss: NDArray[np.float64],
    progress: bool = True,
    **kw,
) -> list[dict]:
    """Cartesian (alpha, S) grid. Returns a flat list of rows."""
    rows = []
    total = len(alphas) * len(Ss)
    for i, a in enumerate(alphas):
        for j, s in enumerate(Ss):
            rows.extend(condition(float(a), float(s), **kw))
            if progress:
                print(f"  [{i * len(Ss) + j + 1}/{total}] alpha={a:.3g} S={s:.3g}", flush=True)
    return rows
