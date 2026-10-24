# TraTT — Physics conventions

Conventions, assumptions and analytic limits behind the solver. The tests in
`tests/` assert the analytic facts below.

Design rule: the solver produces `I(η,t)` and nothing else; every result is an
analysis of that cube.

## 0. Physics conventions

Two-beam Takagi–Taupin, one traveling acoustic modulation, stroboscopic
(quasi-static) X-ray response.

Justification for stroboscopic: the retardation group is `β = Ω·L/v_g ≈ Ω·L/c`.
At f = 1 GHz, L = 10 µm → β ≈ 2×10⁻⁴. Measured once in `tests/`, then dropped.
If a later result needs it, the channel-coupled form is a known extension.

### Fields

Normalized depth `ζ = z/L ∈ [0,1]`, slow time `t̃ = t/t_tr`, `t_tr = L/v_s`.

```
dD0/dζ = iπα · e^{+iφ(ζ;t̃)} · Dh
dDh/dζ = iπα · e^{-iφ(ζ;t̃)} · D0 + iη · Dh
```

### Drive

```
φ(ζ; t̃) = S · A(t̃ − ζ) · sin(q·ζ − w·t̃)
```

`A` = causal envelope, `A(s<0) = 0` gives the wavefront cutoff `z_max = v_s t`.

### Assumptions baked into "two free parameters"

`chi_0` is not assumed equal to anything — it is gauged away exactly. Substitute
`D_0 = e^{kz} d_0` and `D_h = e^{kz} d_h` with `k = i pi chi_0 / (lambda gamma)`,
the *same* factor for both beams; the `gamma k` term then cancels the `chi_0` term
in each equation. This works only because `chi_0` enters both equations with the
same coefficient, i.e. only in symmetric geometry (`gamma_0 = gamma_h`).

Be precise about complex `chi_0`: `Re chi_0` is a pure phase and does vanish from
`|D_h|^2`, but `Im chi_0` leaves a common `exp(-mu z)` that survives in the
absolute intensity. It cancels only in normalized observables — `p_n`, curve
shape, formation times — which is all we use. (This model sets `chi_0 = 0`
outright, so the point is moot here and matters only when comparing to data.)

`chi_0` does not vanish physically. It shifts the Bragg angle by
`|Re chi_0| / sin(2 theta_B)` — about 6.5 arcsec for Si(111) at 8 keV. That shift
is absorbed into the definition of `eta`, so `eta = 0` is the *refraction
corrected* Bragg angle, an offset to remember when matching real rocking data.

The real hidden assumption is in the coupling. It is properly
`C sqrt(chi_h chi_hbar)`, and writing it as `C|chi_h|` requires
`chi_hbar = chi_h*`, i.e. a **non-absorbing** crystal. That is exactly what
`alpha: float` in `core.py` encodes; with absorption `alpha` would be complex.

Full list:

1. two-beam
2. **symmetric geometry** `gamma_0 = gamma_h` — the `chi_0` gauge depends on it
3. plane wave, monochromatic
4. **no absorption**, `Im chi = 0`, so `alpha` is real
5. `chi_h chi_hbar = |chi_h|^2` (follows from 4; Si is centrosymmetric so the
   structure factor causes no extra trouble)
6. fixed polarization factor `C` (sigma: 1, pi: `cos 2 theta_B`), no mixing
7. the acoustic wave enters **only as a phase** — `chi_h` itself is unchanged
   (no elasto-optic change, no vibrational Debye–Waller reduction)

What 4 discards: uniform absorption (harmless — common factor) and the
**Borrmann effect** (not harmless — `Im chi_h` gives the two Bloch modes
different absorption, and mode-selective damping is itself a candidate extra
timescale). Excluded from the Laue solver. The Bragg solver takes it as
`kappa_h` (`c -> i pi alpha + kappa_h / 2`, centrosymmetric crystal; notebook
Entry 77), next to the normal absorption `kappa`.

Note that strain is *not* discarded by 7: `d phi/d zeta = h . du/d zeta` is the
local lattice-spacing change, so it is already inside `phi`.

### Dimensionless groups — only two are free

| symbol | meaning | value |
|---|---|---|
| `α = L/Λ_ext` | dynamical diffraction strength | swept |
| `S = h·u₀` | acoustic modulation index | swept |
| `q = K·L` | carrier spatial phase across slab | `= 2πα_ac` |
| `w = Ω·t_tr` | acoustic phase per transit | `= q` (non-dispersive) |

`α_ac = L/Λ_ac` is the acoustic thickness parameter. **`q = w` is locked** by
`Ω = K v_s`; treating them as independent is a generalization not needed here.

### Analytic facts used as ground truth (not fitted)

1. Unmodulated symmetric Laue, η=0: `|Dh(1)|² = sin²(πα)`.
2. Lossless Laue: `|D0|² + |Dh|² = 1` at every ζ.
3. Phase grating expansion `e^{iφ} = Σ_n J_n(S) e^{in(qζ − w t̃)}` puts
   **channel n at `η_n = n·q`**. Satellite positions are known, never searched for.
4. Weak drive: `p_{±1}/p_0 → (J_1(S)/J_0(S))² → (S/2)²`.
5. Step envelope, t̃ ≥ 1: `φ` is already in its final form, so `I(η,t̃)` is exactly
   periodic in t̃ with period `2π/w`. **`t_f ≤ t_tr` is therefore structural.**
   Any run reporting `t_f > t_tr` under a step envelope is a bug, not a discovery.

Fact 5 reframes the central question of the abstract:

> not "does the transient outlive the wavefront?"
> but **"what fraction of the transit does formation actually take, and why?"**

### Standing prediction

Bragg geometry reflects only from the top ~`Λ_ext`. The front only has to cross
that depth for the diffracted state to be established:

```
t_f / t_tr  ≈  min(1, 1/α)      (Bragg)
t_f / t_tr  ≈  1                (Laue, whole thickness contributes)
```

Quantified in `explore/p0_equations.py`: for the **unmodulated static** problem,
the depth at which `|R|²` reaches 99 % of its final value is

| α | 0.5 | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|
| ζ(99 %) | 0.983 | 0.866 | 0.480 | 0.244 | 0.127 | 0.060 |
| 1/α | 2.00 | 1.00 | 0.50 | 0.25 | 0.125 | 0.062 |

so `ζ_sat = min(1, 1/α)` holds to 2–4 % for α ≥ 4.

This is a static penetration depth, not a formation time — the real problem is
modulated and time dependent, so it may not carry over. That is exactly why it is
useful: a quantitative baseline fixed *before* the sweep, which the parameter sweep either
confirms or breaks. If it holds, the "additional dynamical timescale" of the
abstract is the extinction depth mapped into time, `t_ext = Λ_ext/v_s = t_tr/α`.
