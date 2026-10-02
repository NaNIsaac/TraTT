# Lab notebook — transient formation of acoustic Bragg satellites

Written as the working record of the paper. Every entry states the hypothesis,
the exact setup (script + parameters), the numbers that came out, and what they
do and do not support. Negative results stay in. Scripts live in `explore/nb/`.

Model throughout unless stated: stroboscopic two-beam Takagi–Taupin, Laue,
non-absorbing, step onset (rise = 0), resolved channels (alpha_ac = 20,
q = 125.66, so 2 pi alpha / q = alpha / 20).

Notation: t = t_real / t_tr (filled depth fraction), alpha_n = alpha |J_n(S)|.

---

## 2026-09-30 — Entry 0: starting point

Established (handoff C-1.10, commit c8c7193): satellite n during formation is a
Laue slab of thickness t and coupling pi alpha_n.

    I(nq, t) ~ sin^2(pi alpha_n t)
    p_n(t)   ~ pi^2 alpha_n int_0^{2 pi alpha_n t} J_0(x) dx

Open questions queued for this session:

1. Formation order. Does the order in which channels reach their first maximum
   follow |J_n(S)|, and can S reverse it?
2. The +12 % excess of p_1 near its first maximum at S = 0.5, alpha = 2.
3. The negative shift of the satellite peak position.
4. Spectral width: does it cross over from the finite-depth 5.566/t to a
   dynamical width set by alpha_n?
5. Finite rise time: what survives of the temporal Pendellösung?
6. Units: is the effect observable at a real source?
7. Prior art for the satellite extinction length Lambda / |J_n(S)|.

---

## 2026-09-30 — Entry 1: formation order and the n = 0 two-layer law

**Hypothesis.** (a) Channel n != 0 reaches its first centre maximum at
t* = 1 / (2 alpha |J_n(S)|), so the order is set by |J_n(S)|, not by |n|.
(b) The n = 0 channel is not isolated: the filled layer (coupling pi alpha J_0)
and the unfilled layer (pi alpha) are both on Bragg at eta = 0. Laue transfer
matrices at exact Bragg are rotations, so the angles add:

    I(0, t) = sin^2( pi alpha [1 - (1 - J_0(S)) t] )

**Setup.** `explore/nb/e1_order.py`. alpha_ac = 20, n_z = 1200, 201 times in
[0, 1], intensity read at eta = n q exactly. (S, alpha) = (0.5, 2), (1, 2),
(3, 2), (3, 1).

**Result.** max |I_sim - I_th| over t in [0, 1]:

| S, alpha | n=0 | n=1 | n=2 | n=3 |
|---|---|---|---|---|
| 0.5, 2 | 0.005 | 0.079 | 0.008 | 0.000 |
| 1, 2   | 0.015 | 0.037 | 0.048 | 0.005 |
| 3, 2   | 0.028 | 0.055 | 0.060 | 0.058 |
| 3, 1   | 0.013 | 0.024 | 0.028 | 0.030 |

First-maximum time where it falls inside the transit (sim vs prediction):

| S, alpha | n | J_n(S) | t* sim | t* pred |
|---|---|---|---|---|
| 1, 2 | 1 | 0.440 | 0.570 | 0.568 |
| 3, 2 | 2 | 0.486 | 0.500 | 0.514 |
| 3, 2 | 1 | 0.339 | 0.750 | 0.737 |
| 3, 2 | 3 | 0.309 | 0.825 | 0.809 |

**Reading.**
- (a) holds to within the time grid plus ~2 %. At S = 3 the channels peak
  in the order n = 2, 1, 3: **formation order is not ordered by harmonic
  index.** It follows the Bessel ladder, so sweeping S reorders it. This is a
  prediction a time-resolved experiment can test channel by channel.
- (b) holds to 0.005–0.028. Two consequences worth stating in the paper:
  the main line intensity oscillates in time with frequency alpha (1 - J_0(S))
  during transit, and for J_0(S) < 0 (S = 3) the filled layer *counter-rotates*,
  so the main line can pass through zero while the front is inside the crystal.
- Largest n = 1 deviation (0.079 at S = 0.5, alpha = 2) is the same case as the
  +12 % of C-1.10 — see Entry 2.

**Does not show.** Nothing about the *integrated* order; only centre values.

---

## 2026-09-30 — Entry 2: the +12 % excess is the static rocking-curve tail

**Hypothesis.** The excess of p_1 over the collapse law is the unmodulated
crystal's off-Bragg tail landing in the n = 1 window. Off Bragg,
|Dh|^2 ~ (2 pi alpha / eta)^2 sin^2(eta / 2); fringe-averaged over the window
(q/2, 3q/2) this integrates to

    p_tail ~ 8 (pi alpha)^2 / (3 q)

For S = 0.5, alpha = 2, alpha_ac = 20: p_tail = 0.84, p_1 at first maximum ~ 7.0,
ratio 12 % — the size of the excess seen in C-1.10.

**Setup.** `explore/nb/e2_tail.py`. S = 0.5, n = 1, alpha in {1, 2},
alpha_ac in {10, 20, 40}, t in {0.25, 0.5, 0.75, 1.0}. Window integral of raw I
minus the collapse law = "excess". Compare with the estimate and with the
window integral of the S = 0 crystal.

**Result.**

| alpha | alpha_ac | estimate | S=0 window | excess (t = .25/.5/.75/1) | residual after S=0 subtraction |
|---|---|---|---|---|---|
| 1 | 20 | 0.209 | 0.210 | 0.222 0.220 0.216 0.260 | 1.7 % |
| 1 | 40 | 0.105 | 0.105 | 0.110 0.109 0.108 0.131 | 0.9 % |
| 2 | 20 | 0.838 | 0.839 | 0.866 0.848 0.807 0.962 | 1.9 % |
| 2 | 40 | 0.419 | 0.419 | 0.437 0.430 0.420 0.509 | 1.4 % |
| 1 | 10 | 0.419 | 0.421 | 0.208 0.434 0.176 0.498 | 24 % |
| 2 | 10 | 1.676 | 1.604 | 0.733 1.504 0.615 1.574 | 25 % |

**Reading.**
- Confirmed for alpha_ac >= 20: the excess is the static tail, scales as
  alpha^2 / alpha_ac, and the closed estimate matches the S = 0 control to
  < 0.5 %. Subtracting the S = 0 window leaves the collapse law good to 1–2 %.
  The collapse figure should plot p_n - p_n^{S=0}.
- At alpha_ac = 10 the excess alternates with t (low at t = .25, .75; high at
  .5, 1). The tail is not fringe-averaged there and interferes coherently with
  the satellite. **Open.** Period 0.5 in t at alpha_ac = 10 — check whether it is
  5 / alpha_ac or tied to q in some other way before interpreting.
- At t = 1 the residual is consistently larger (+0.03 to +0.12). The unfilled
  layer is gone there, so the tail comes from the filled crystal (coupling
  alpha J_0 in the n = 0 channel) plus weak neighbours. Not resolved; minor.

**Follow-up (same day).** Fine time grid (361 points), S = 0.5, alpha = 1,
excess minus S = 0 window, FFT in t:

| alpha_ac | dominant frequency (per unit t) | rms amplitude |
|---|---|---|
| 7  | 6.65 | 0.236 |
| 10 | 9.97 | 0.179 |
| 14 | 14.4 | 0.127 |

The alternation is a ripple at the **acoustic period 1/alpha_ac**, with
amplitude falling roughly as 1/alpha_ac: satellite–tail interference that
follows the drive phase q t (micromotion). It was invisible at alpha_ac = 20 and
40 only because t = .25/.5/.75/1 are whole numbers of acoustic periods there —
an aliasing accident, not an absence. **Lesson for all later entries: sample t
off the acoustic period grid, or average over one period, before reading a
slow envelope.** The collapse law describes the period-averaged envelope.

---

## 2026-09-30 — Entry 3: the satellite peak shift is kinematic interference (negative result)

**Question.** C-1.10 noted a negative shift of the n = 1 peak (up to -5 at
alpha = 4). Is it a dynamical resonance shift (Floquet level shift,
~ (pi alpha J)^2 / q) or an artefact of reading a peak on top of other
contributions?

**Setup.** `explore/nb/e3_shift.py` (t = 0.3, 0.6, 1.0) and
`explore/nb/e3b_shift_scan.py` (t = 0.3137, chosen off the acoustic period grid).
Peak position = parabola through the three highest samples, eta step 0.01.
Both raw |A|^2 and the amplitude difference |A - A0|^2. S = ±1.5,
alpha_ac in {20, 40, 80}, alpha in {0.01, 0.1, 0.5, 1}, n = ±1.

**Result (t = 0.3137, raw, n = +1 / -1).**

| S | alpha_ac | alpha = 0.01 | alpha = 0.5 | alpha = 1 |
|---|---|---|---|---|
| +1.5 | 20 | +0.25 / +2.82 | +0.22 / +2.83 | +0.14 / +2.85 |
| +1.5 | 40 | +1.47 / +1.23 | +1.44 / +1.23 | +1.35 / +1.25 |
| +1.5 | 80 | -0.57 / -0.95 | -0.56 / -0.93 | -0.55 / -0.84 |
| -1.5 | 20 | -2.82 / -0.25 | -2.83 / -0.22 | -2.85 / -0.14 |

**Reading.**
- The shift is the same at alpha = 0.01 as at alpha = 1. It is **kinematic**, not
  dynamical. The Floquet estimate (pi alpha J_0)^2 / q is 0.02 at alpha = 1 —
  two orders smaller than what is measured.
- Its sign and size jump with alpha_ac at fixed t: the shift follows a phase,
  not a coupling strength. Consistent with interference between the satellite
  amplitude and the unmodulated off-Bragg amplitude A0 near eta = n q. In the
  kinematic limit both scale as alpha, so their ratio — and the shift — is
  alpha-independent, exactly as observed.
- S -> -S maps (n, shift) -> (-n, -shift), the symmetry of the drive.
- The shift decays toward t = 1 (Entry 3 first table: -1.6, -0.5, -0.1 at
  t = 0.3, 0.6, 1 for alpha = 0.5), as the satellite grows relative to A0.
- The large values at alpha = 2, t = 1 (-5.9) are a different artefact:
  sin^2(pi * 1.12) is near zero, the centre dips and the parabola lands on a
  side lobe.

**Conclusion.** Peak position is **not** a dynamical observable in this model.
It must not go on the poster as a result. If the abstract's "peak position" is
kept, report it as "fixed at n q to within the finite-depth width; apparent
shifts of order 1–3 in eta are kinematic interference with the unmodulated
crystal". Removing A0 (difference observable) does not cure it either: a
residual +0.2 to +0.4 at alpha -> 0 remains, from neighbour-channel sinc tails.

---

## 2026-09-30 — Entry 4: satellite width narrows below the finite-depth limit

**Hypothesis.** If satellite n is a Laue slab (thickness t, coupling
a = pi alpha_n), its line shape is

    f(d) = a^2 sin^2( t sqrt(a^2 + d^2/4) ) / (a^2 + d^2/4),   d = eta - n q

so t * FWHM = G(alpha_n t) is a universal function. G(0) = 5.566 is the
finite-depth sinc^2 of C-1.8 ("width = 1/t"). The dynamical question is how G
departs from 5.566.

**Setup.** `explore/nb/e4_width.py`. G computed with the same `peak_fwhm`
routine used on the data. alpha_ac = 20, n_z = 1200, eta step 0.03,
46 times in [0.1, 1] offset by 0.0037 (off the acoustic grid). Raw I, and the
amplitude difference |A - A0|^2 (`satellite_cube`) as a second observable.

**Theory curve.**

| alpha_n t | 0.02 | 0.30 | 0.58 | 0.72 | 0.86 | 1.00 | 1.14 | 1.28 |
|---|---|---|---|---|---|---|---|---|
| t * FWHM | 5.57 | 5.38 | 4.79 | 4.23 | 3.96 | 4.48 | 5.61 | 13.96 |

**Result.**
- Raw I, alpha_n t >= 0.5: matches G to ~1–2 % in all four cases,
  (S, alpha, n) = (1.5, 1, 1), (1.5, 2, 1), (3, 2, 2), (0.5, 4, 1).
  Example (1.5, 2, 1): alpha_n t = 0.65 / 0.79 / 0.92 / 1.05 gives
  4.48 / 3.95 / 4.11 / 4.74 against 4.54 / 3.85 / 4.14 / 4.79.
- Raw I, alpha_n t < 0.4: off by up to 50 % (low). This is the C-1.8 problem:
  the unfilled crystal's full-thickness fringes sit inside the satellite.
- Difference observable, early times: 5.53 / 5.43 / 5.37 against 5.56 / 5.55
  (case 1.5, 1, 1), good to 1–3 %, and max deviation 3.8 % for alpha_n t < 0.9
  in that case. For the stronger cases (alpha = 2, 4) its worst point is 18 % /
  33 %, located later in the transit — the raw observable is the better one there.

**Reading.**
- The spectral width has a **dynamical minimum**: t * FWHM falls to 3.96, i.e.
  0.71 of the finite-depth value, at alpha_n t ~ 0.86, then broadens again and
  diverges as the centre drops to zero (peak splitting) at alpha_n t = 1.
- So the width is *not* transit-limited 1/t throughout. In time,
  W(t) = G(alpha_n t) / t narrows faster than 1/t between alpha_n t ~ 0.4 and
  0.86, then stops narrowing — an observable signature of the dynamical
  timescale t_tr / alpha_n that the integrated intensity only shows as a hump.
- This resolves C-1.8 in part: width is measurable, with the observable chosen
  by regime (difference early, raw late) and compared to G, not to 5.566/t.

**Does not show.** Behaviour past the first zero (alpha_n t > 1) was not
examined; there the line is double-peaked and FWHM is the wrong statistic.

---

## 2026-09-30 — Entry 5: finite rise time — the "accumulated coupling" law

**Hypothesis.** With a rise time, depth zeta carries modulation S A(t - zeta).
If the envelope is slow on the acoustic wavelength, harmonic n has local
coupling pi alpha J_n(S A(t - zeta)). At exact resonance the Laue transfer
matrices of successive layers are rotations about the same axis, so the angles
add:

    I(n q, t) = sin^2( pi alpha Theta_n(t) ),   Theta_n(t) = | int_0^1 J_n(S A(t - zeta)) dzeta |

This one formula contains everything so far: rise = 0 gives Theta_n = |J_n| t
(n != 0, C-1.10) and Theta_0 = 1 - (1 - J_0) t (Entry 1).

**Setup.** `explore/nb/e5_rise.py`. alpha_ac = 20, n_z = 1500, 121 times in
[0, 1 + rise] offset by 0.0037, centre intensity at eta = n q, n = 0, 1, 2.
Comparison: the same formula with a step envelope ("step law").

**Result.** max |I_sim - law| over the whole window:

| S, alpha, rise | n=0 | n=1 | n=2 | step law (worst n) |
|---|---|---|---|---|
| 1.5, 2, 0.05 | 0.024 | 0.035 | 0.050 | 0.103 |
| 1.5, 2, 0.3  | 0.023 | 0.036 | 0.054 | 0.535 |
| 1.5, 2, 0.8  | 0.024 | 0.042 | 0.051 | 0.994 |
| 3, 2, 0.5    | 0.022 | 0.052 | 0.058 | 0.996 |
| 1.5, 4, 0.5  | 0.129 | 0.181 | 0.191 | 0.988 |

alpha = 4 case versus alpha_ac (n_z = 60 alpha_ac):

| alpha_ac | n=0 | n=1 | n=2 |
|---|---|---|---|
| 20 | 0.129 | 0.181 | 0.191 |
| 40 | 0.020 | 0.061 | 0.045 |
| 80 | 0.010 | 0.013 | 0.020 |

**Reading.**
- The accumulated-coupling law holds for every rise tried, including
  rise = 0.05 = one acoustic period, where "slow envelope" is not obviously
  satisfied. The step law fails completely once rise >~ 0.3.
- The residual falls roughly as (alpha / alpha_ac)^2: the law is the leading
  term of an expansion in pi alpha / q, and the corrections are the
  off-resonant channels.
- For S = 3 the n = 0 coupling J_0(S A) changes sign inside the ramp; the law
  still holds because the signed integral is taken before the modulus.
- **Physical meaning for the paper:** a finite rise does not add a new
  timescale of its own; it re-weights the depth integral. The formation curve
  is sin^2 of the *accumulated* coupling seen along the depth, and the rise only
  changes how fast that accumulation proceeds. This is also why C-1.7(b) found
  the overshoot ratio insensitive to rise: the maximum of sin^2 is 1 whatever
  path Theta takes to pi/2, only the time of the peak moves.

---

## 2026-09-30 — Entry 6: units — the temporal Pendellösung period is thickness-independent

**Observation (algebra, then numbers).** In real time the satellite period is

    tau_n = t_tr / alpha_n = (L / v_s) / (alpha |J_n(S)|) = Lambda_P / (v_s |J_n(S)|)

with Lambda_P = L / alpha the intensity Pendellösung period (the code's
convention, handoff 7.3). **L cancels.** tau_n is the time the acoustic front
needs to cross one *satellite* Pendellösung length Lambda_P / |J_n(S)|. The
crystal thickness only sets how many oscillations fit in the transit
(alpha_n of them).

**Numbers.** `explore/nb/e6_units.py`. Inputs are textbook-level and are
flagged in the script as TO VERIFY with X0h/XOP before publication:
Lambda_P(Si 220, ~8 keV, sigma) ~ 15 um; shear wave in Si [100]
v_T = 5.84 km/s (a shear wave is needed so that h.u != 0 in symmetric Laue,
handoff 7.4); d_220 = 1.920 A. alpha_ac chosen >= max(20, 5 alpha) to keep the
channels resolved.

| L | alpha | t_tr | f_ac | S | u0 | strain | alpha_1 | tau_1 | first max |
|---|---|---|---|---|---|---|---|---|---|
| 50 um  | 3.3  | 8.6 ns  | 2.3 GHz | 0.5 | 15 pm | 3.8e-5 | 0.81 | 10.6 ns | 5.3 ns |
| 50 um  | 3.3  | 8.6 ns  | 2.3 GHz | 1.5 | 46 pm | 1.2e-4 | 1.86 | 4.6 ns | 2.3 ns |
| 100 um | 6.7  | 17 ns   | 2.0 GHz | 1.5 | 46 pm | 9.6e-5 | 3.72 | 4.6 ns | 2.3 ns |
| 200 um | 13.3 | 34 ns   | 2.0 GHz | 1.5 | 46 pm | 9.6e-5 | 7.44 | 4.6 ns | 2.3 ns |

**Reading.**
- Timescales of a few ns, against ~50–100 ps synchrotron bunches: the
  oscillation is resolvable at a storage ring with a gated GHz transducer
  synchronised to the bunch clock. No XFEL needed.
- Strain ~1e-4 at S = 1.5 is within reach of GHz piezo transducers (order of
  magnitude only; the achievable value is transducer-specific).
- Because tau_n = Lambda_P / (v_s |J_n(S)|), measuring tau_1 and tau_2 at one
  drive setting gives |J_2(S) / J_1(S)| — **an in-situ calibration of the
  acoustic amplitude S that does not need the absolute intensity**. Worth one
  sentence in the paper; a ratio of periods is robust against detector gain,
  flux and absorption.
- GHz with alpha_ac >= 20 means Lambda_ac of 2.5–5 um. Whether such a shear
  wave can be launched cleanly through a 50–200 um slab is an experimental
  question this model does not address (acoustic attenuation and diffraction
  are absent).

**Caveats.** Absorption is absent (mu L for Si 220 at 8 keV and L = 100 um is
~1.4, not negligible; the Borrmann effect makes it mode-selective). The numbers
above are for planning, not for a claim.

---

## 2026-09-30 — Entry 7: prior art — what is already known, what is left

Web search only (abstracts and summaries, full texts NOT read). Every item
below must be checked in the original before it is cited.

**Already known — must be cited, cannot be claimed.**
- **Temporal Pendellösung of the main reflection.** DeCamp, Reis, Bucksbaum
  et al., "Coherent control of pulsed X-ray beams", Nature 413, 825 (2001):
  a laser-generated strain pulse entering a Laue crystal (Borrmann switch)
  is described as "the temporal analog of the Pendellösung effect", the pulse
  acting "like a thin interface layer dividing the crystal"; switching needs
  the pulse to travel ~1/4 Pendellösung depth. Follow-up: DeCamp et al.,
  PRL 91, 165502 (2003) (Ge, anomalous transmission). A summary of this line of
  work states that diffracted intensity oscillates in time with period
  Pendellösung length / sound speed.
  → Our n = 0 two-layer law (Entry 1) is the same physics for a CW grating.
  Their dividing layer is a pulse; ours is the front of a filled grating whose
  n = 0 coupling is J_0(S). The rotation-angle-addition statement and the
  counter-rotation for J_0 < 0 are the CW generalisation, not a new effect.
- **Satellite extinction length in the steady state.** Work on X-ray
  diffraction by ultrasonic superlattices (Entin and co-workers; Zolotoyabko
  and co-workers, 1970s–2000s) states that ultrasound creates satellite
  reflections with a very large effective extinction length. The J_n(S)
  scaling of that length is standard two-beam Floquet algebra and is almost
  certainly in that literature.

**Not found in this search — candidate contributions, to be defended.**
1. The **time-domain mapping** for a CW travelling grating: satellite n
   during transit equals a Laue slab of thickness v_s t and coupling
   pi alpha |J_n(S)|, verified against full TT with **no fit parameter**, for
   centre intensity, integrated intensity and width (C-1.10, Entries 2, 4).
2. The **collapse** of all (alpha, S, n) onto one curve of alpha_n t.
3. **Formation order set by the Bessel ladder**, reorderable with S
   (Entry 1: n = 2 before n = 1 before n = 3 at S = 3).
4. **Accumulated-coupling law** for finite rise time (Entry 5): one formula,
   sin^2(pi alpha | int J_n(S A(t - zeta)) dzeta |), covering step and ramp
   onsets, all channels.
5. **Dynamical narrowing** of the satellite width to 0.71 of the
   finite-depth limit at alpha_n t ~ 0.86 (Entry 4).
6. **In-situ calibration of S** from the ratio of temporal periods,
   tau_2 / tau_1 = |J_1(S) / J_2(S)|, thickness- and flux-independent (Entry 6).

**Framing for the poster.** Do not sell "temporal Pendellösung" as new; sell
it as the organising principle that turns the transient satellite problem
into a solved one — and list 2–6 as what follows.

Sources: nature.com/articles/35101560 ; arxiv.org/abs/cond-mat/0301002 ;
osti.gov/etdeweb/biblio/20231834 (ultrasonic superlattice, full external
reflection) ; arxiv.org/abs/2302.14116 (picosecond ultrasonics with x-rays).

**Entry 2b (same day) — the tail correction does not generalise.** Residual of
the integrated law in units of pi^2 alpha_n, 100 times off the period grid:

| S, alpha, n | raw rms / max | minus S=0 window | minus S=0 * Theta_0(t)^2 |
|---|---|---|---|
| 0.5, 2, 1 | 0.145 / 0.264 | 0.076 / 0.161 | 0.070 / 0.149 |
| 1.5, 2, 1 | 0.044 / 0.105 | 0.106 / 0.181 | 0.079 / 0.156 |
| 1.5, 4, 1 | 0.081 / 0.185 | 0.202 / 0.323 | 0.152 / 0.290 |
| 2, 4, 2   | 0.068 / 0.219 | 0.052 / 0.172 | 0.058 / 0.174 |
| 1, 6, 1   | 0.085 / 0.302 | 0.257 / 0.462 | 0.209 / 0.422 |

Subtracting the static tail helps only when alpha_n is small relative to the
tail (first row). For strong satellites it over-subtracts: during transit the
tail is not the static one (the n = 0 coupling in the filled layer is
alpha J_0, and the tail interferes with the satellite rather than adding).
Scaling the tail by Theta_0^2 goes the right way but not far enough. **Decision:**
the collapse figure uses raw window integrals; the integrated law is quoted as
good to ~5–15 % of pi^2 alpha_n (rms), the centre law as the precise one.
The Entry 2 conclusion stands only for the case it was measured in
(small alpha_n, alpha_ac >= 20).

---

## 2026-09-30 — Entry 8: the transient ±n asymmetry lives in the first acoustic wavelength

**Question.** C-1.7(c) reported a transient asymmetry
A_1 = (p_+1 - p_-1)/(p_+1 + p_-1) with max |A| ~ 0.5, independent of alpha, and
called it a geometric effect of the travelling front. The collapse law is
symmetric in n. Is the asymmetry robust as channels become better resolved?

**Setup.** `explore/nb/e8_asym.py`. S = 0.5, alpha in {0.05, 1.3},
alpha_ac in {5, 10, 20, 40}, 199 times in (0, 1), window integrals of raw I.

**Result.**

| alpha | alpha_ac | max abs A | t at max | t * alpha_ac | intensity-weighted rms | rms for t > 0.3 |
|---|---|---|---|---|---|---|
| 0.05 | 5  | 0.454 | 0.160 | 0.80 | 0.156 | 0.126 |
| 0.05 | 10 | 0.510 | 0.080 | 0.80 | 0.104 | 0.080 |
| 0.05 | 20 | 0.501 | 0.040 | 0.80 | 0.060 | 0.043 |
| 0.05 | 40 | 0.496 | 0.020 | 0.80 | 0.034 | 0.022 |
| 1.3  | 5  | 0.507 | 0.150 | 0.75 | 0.161 | 0.119 |
| 1.3  | 10 | 0.503 | 0.080 | 0.80 | 0.104 | 0.074 |
| 1.3  | 20 | 0.499 | 0.040 | 0.80 | 0.062 | 0.040 |
| 1.3  | 40 | 0.496 | 0.020 | 0.80 | 0.035 | 0.021 |

**Reading.**
- The peak asymmetry is always ~0.5, and it always occurs at
  t * alpha_ac ~ 0.8: when the filled depth holds **less than one acoustic
  wavelength**. At that point the "channels" are not yet defined — the
  finite-depth width 5.566/t is larger than the channel spacing q.
- The weighted rms falls as 1/alpha_ac (0.156, 0.104, 0.060, 0.034 on
  doubling alpha_ac), and is alpha-independent: kinematic, like Entry 3.
- **Correction to C-1.7(c).** The asymmetry is real in the model but its clock
  is the number of filled acoustic wavelengths, alpha_ac t, not the transit
  and not dynamics. Its "4–6x larger than steady state" statement is a
  property of alpha_ac = 5. It should not be presented as a dynamical feature.
  It can be presented, if at all, as: "for the first acoustic wavelength the
  satellites are not yet symmetric sidebands; symmetry is restored on a time
  1/alpha_ac (one acoustic period)".

---

## 2026-09-30 — Entry 9: Bragg geometry — the original t_f ~ 1/alpha hypothesis, corrected to 1/alpha_n

**Hypothesis.** Same mapping, Bragg slab: the front enters at the entrance
face, the filled layer [0, t] is a Bragg slab of coupling pi alpha_n at
eta = n q, the unfilled remainder is off-Bragg. Non-absorbing Bragg slab at
exact Bragg: I(n q, t) = tanh^2(pi alpha_n t). No oscillation — saturation.

**Setup.** `explore/nb/e9_bragg.py`. alpha_ac = 20, n_z = 1500,
geometry = "bragg", centre intensity. Then 2001 times to locate the 99 % point
of the final value.

**Result.** max |I - tanh^2(pi alpha_n t)| over [0, 1]:
0.004 (S 1.5, alpha 0.5, n 1), 0.006 (1.5, 1, 1), 0.031 (1.5, 2, 1),
0.043 (3, 2, 2).

99 % formation time, prediction t99 = atanh(sqrt(0.99) tanh(pi alpha_n)) / (pi alpha_n)
(handoff 7.3 formula with alpha -> alpha_n):

| S, alpha, n | alpha_n | t99 sim | t99 pred |
|---|---|---|---|
| 1.5, 1, 1 | 0.56 | 0.974 | 0.977 |
| 1.5, 2, 1 | 1.12 | 0.813 | 0.810 |
| 3, 2, 2   | 0.97 | 0.863 | 0.876 |
| 1.5, 4, 1 | 2.23 | 0.413 | 0.427 |

**Reading.**
- The mapping holds in Bragg too. The satellite forms by *saturation*
  (tanh^2) rather than oscillation (sin^2); which one you get is the geometry.
- **The handoff's original hypothesis t_f / t_tr ~ min(1, 1/alpha) (section 9.3)
  is right in form and wrong in the variable.** The satellite's own extinction
  sets the time: t99 ~ 0.953 / (alpha |J_n(S)|) for large alpha_n, capped at 1.
  With alpha_n rather than alpha, it holds to ~1–3 %.
- Why session 1 saw no 1/alpha: in Laue there is no saturation to reach early,
  and the formation-time definition keyed on the full-profile residual always
  ends at t = 1. In Bragg the centre intensity genuinely reaches its final
  value before the transit ends.
- So the answer to the abstract's central question differs by geometry:
  **Bragg — the satellite formation time is dynamical, 0.95 t_tr / alpha_n;
  Laue — the transient is a dynamical oscillation of period t_tr / alpha_n and
  only the end point is set by the transit.**

**Does not show.** Only centre intensity. Bragg integrated intensity and width
not checked. Bragg is deprioritised by the user — this was a cheap check to
test generality, not a full study.

---

## 2026-09-30 — Entry 10: poster summary figure

`explore/make_poster_figure.py` -> `figures/poster_mapping.png`. Four panels,
closed-form curves, full-TT points, alpha_ac = 20, times off the period grid:

- a Laue centre intensity of 4 (S, alpha, n) cases on sin^2(pi alpha_n t).
- b Bragg, same cases, on tanh^2(pi alpha_n t).
- c First-maximum time vs S for n = 1, 2, 3 at alpha = 2, theory
  1 / (2 alpha |J_n(S)|) and simulated points; S = 3 shows the n = 2, 1, 3 order.
- d Rise = 0.5, S = 1.5, alpha = 2, n = 0, 1, 2 on the accumulated-coupling law,
  with the step-onset law for n = 1 dashed — visibly wrong.

Palette: reference categorical order (blue, orange, aqua, yellow), validated
with the dataviz validator (all PASS; aqua contrast WARN, relieved by marker
shape + legend on every series).

Worth noticing in panel d: after t = 1 + rise the channels sit at their steady
values sin^2(pi alpha |J_n(S)|) — here n = 2 near 1 and n = 0 near 0 — so the
*steady* ladder is the end point of the same law, not a separate result.

---

## 2026-09-30 — Entry 11: absorption — estimate only, NOT simulated

The code is non-absorbing (real alpha). Estimate for planning:

- A common chi_0'' gives exp(-mu L) overall and cancels in shapes.
- The Borrmann part (chi_h'') gives the two Laue modes absorption
  mu (1 ± kappa); in the filled layer the satellite's coupling is scaled by
  J_n, so its modes split by mu kappa |J_n(S)|. The sin^2 then becomes
  [cosh(mu L kappa |J_n| t) - cos(2 pi alpha_n t)] / 2 up to an overall
  exp(-mu L) factor, and the oscillation contrast is ~ 1 / cosh(mu L kappa |J_n| t).
- Si 220, ~8 keV, L = 100 um: mu L ~ 1.4, kappa ~ 0.95 (TO VERIFY), J_1(1.5) = 0.56,
  so at t = 1 the argument is ~0.75 and the contrast ~0.78. **The oscillation
  survives absorption at the level relevant for the poster's example.**
- Thick crystals (mu L >> 1) kill the contrast; the estimate says keep
  mu L kappa |J_n| <~ 1.

To make this a result: make alpha complex in `core.solve_laue` (the RK4 already
carries complex state; only the alpha >= 0 check and the signature change, as
the ponytail note in core.py says) and redo Entry 1 with a complex coupling.

**Entry 11b (same day) — absorption simulated.** `core.solve_laue` now accepts
a complex alpha (one-line change: the non-negativity check uses the real part;
docstring states the centrosymmetric assumption chi_h = chi_hbar). Tests: 37 pass.
`explore/nb/e11_absorb.py`, centre intensity, alpha_ac = 20, 101 times.

Prediction: the mapping is unchanged with a complex coupling,
I(n q, t) = |sin(pi alpha J_n t)|^2 = [cosh(2 pi alpha'' |J_n| t) - cos(2 pi alpha' |J_n| t)] / 2.

| S, alpha, n | max abs(I - th) | peak sim / th | first minimum sim / th |
|---|---|---|---|
| 1.5, 2+0.05i, 1 | 0.048 | 0.983 / 1.001 | 0.005 / 0.006 |
| 1.5, 2+0.15i, 1 | 0.079 | 0.970 / 1.014 | 0.051 / 0.056 |
| 3, 2+0.10i, 2   | 0.059 | 0.991 / 1.006 | 0.025 / 0.029 |
| 1.5, 4+0.10i, 1 | 0.213 | 0.985 / 1.014 | 0.019 / 0.025 |

The mapping survives Borrmann absorption with the same accuracy as without it
(the alpha = 4 residual is the known alpha/alpha_ac correction of Entry 5).
Link to the estimate above: 2 pi alpha'' = mu L kappa, so for the Si 220,
L = 100 um example alpha'' ~ 0.21 — the range simulated here. The oscillation
minima lift from 0 to ~0.05 at alpha'' = 0.15: contrast stays high.
Still ignored: the common exp(-mu L), and any difference between chi_h and
chi_hbar (non-centrosymmetric crystals).

---

## 2026-09-30 — Entry 12: claims and evidence (draft for the poster / paper)

Status key: **V** verified against full TT here, no fit; **E** estimate only;
**X** withdrawn.

| # | Claim | Evidence | Status |
|---|---|---|---|
| 1 | During transit, satellite n is a crystal of the filled depth t with coupling pi alpha abs(J_n(S)) | C-1.10; Entries 1, 4, 9; `figures/poster_mapping.png` | V (resolved channels, alpha/alpha_ac <~ 0.2) |
| 2 | Laue: centre intensity sin^2(pi alpha_n t) — temporal Pendellösung with period t_tr / alpha_n = Lambda_P / (v_s abs(J_n)) | C-1.10, Entry 1, Entry 6 | V; physics of the main-line analogue is prior art (DeCamp 2001) |
| 3 | Laue: integrated intensity pi^2 alpha_n int_0^{2 pi alpha_n t} J_0 | C-1.10, Entry 2b | V to ~5–15 % rms of pi^2 alpha_n; static tail limits it |
| 4 | Bragg: centre intensity tanh^2(pi alpha_n t); 99 % formation time 0.953 t_tr / alpha_n | Entry 9 | V (centre only) |
| 5 | Formation order set by abs(J_n(S)), reorderable with S (n = 2, 1, 3 at S = 3) | Entry 1, poster panel c | V |
| 6 | n = 0 main line: sin^2(pi alpha [1 - (1 - J_0) t]); counter-rotation when J_0 < 0 | Entry 1 | V; CW analogue of DeCamp's dividing layer |
| 7 | Finite rise: sin^2(pi alpha abs(int_0^1 J_n(S A(t - zeta)) dzeta)) for all n | Entry 5, poster panel d | V, residual ~ (alpha/alpha_ac)^2 |
| 8 | Satellite width narrows to 0.71 of the finite-depth limit at alpha_n t ~ 0.86 | Entry 4 | V for alpha_n t >= 0.5 (raw) and early times (difference observable) |
| 9 | Survives Borrmann absorption at mu L ~ 1.4 | Entry 11b | V (complex alpha); common exp(-mu L) ignored |
| 10 | S calibrated in situ from tau_2 / tau_1 = abs(J_1 / J_2) | Entry 6 | E (algebra; no simulated "measurement" yet) |
| 11 | Timescales ns, strain ~1e-4, GHz: storage-ring accessible | Entry 6 | E (material constants TO VERIFY) |
| — | Peak position shifts dynamically | Entry 3 | X — kinematic interference |
| — | ±n asymmetry is a transit/dynamical effect | Entry 8 | X — first-acoustic-wavelength effect, ∝ 1/alpha_ac |
| — | "No additional timescale" (session 1, C-1.9) | C-1.10 | X — t_tr / alpha_n exists |

**One-sentence story.** A travelling acoustic wave does not switch satellites
on; it grows, behind its front, a crystal whose satellite reflections have
extinction lengths Lambda / abs(J_n(S)) — so each satellite replays its
static thickness curve in time, Pendellösung in Laue and extinction
saturation in Bragg, and the order, timing and width of their appearance are
fixed by the Bessel ladder.

**Open before submission.**
- Read DeCamp 2001 and the ultrasonic-superlattice papers in full; decide the
  exact novelty sentence.
- Verify Lambda_P, kappa and v_T for the chosen reflection with X0h/XOP.
- Width beyond alpha_n t = 1 (split line), and integrated width observable.
- Non-resonant analytic line shape for the whole I(eta, t) near each satellite
  (Entry 4 did the width only).

---

## 2026-09-30 — Entry 13: in-situ S calibration from two periods (claim 10 tested)

**Setup.** `explore/nb/e13_calib.py`. Simulated measurement: centre intensity
of n = 1 and n = 2, 81 times in [0, 1], multiplied by an unknown gain (0.37)
and 3 % multiplicative noise. Each trace fitted with c sin^2(pi f t),
c and f free, multi-start in f. Then S from J_2(S)/J_1(S) = f_2/f_1 on the
branch S < 3.83. **alpha and the gain are not used.**

**Result.**

| S true | alpha | f_1 fit (theory) | f_2 fit (theory) | S recovered | error |
|---|---|---|---|---|---|
| 1.0 | 3 | 1.329 (1.320) | 0.383 (0.345) | 1.094 | +9.4 % |
| 1.5 | 3 | 1.676 (1.674) | 0.707 (0.696) | 1.518 | +1.2 % |
| 2.0 | 3 | 1.726 (1.730) | 1.070 (1.059) | 2.018 | +0.9 % |
| 2.5 | 4 | 1.971 (1.988) | 1.796 (1.784) | 2.519 | +0.8 % |

**Reading.**
- Claim 10 upgraded from E to **V (in simulation)**: ~1 % on S from periods
  alone, gain-free and thickness-free, whenever both channels complete at least
  half an oscillation in the transit (alpha_n >~ 0.5).
- The S = 1 failure is understood: alpha_2 = 0.35, the n = 2 trace shows
  only the rising third of one sin^2, and frequency is degenerate with the
  unknown gain. Rule for an experiment: choose L so that alpha |J_2(S)| >= 0.5.
- Small systematic bias (+1 %) comes from the alpha/alpha_ac corrections
  (fitted f_2 slightly above theory in every row).

---

## 2026-09-30 — Entry 14: a closed three-mode, two-layer model of I(eta, t)

**Model.** `explore/nb/e14_threemode.py`. Near channel n keep, in the filled
layer, the forward beam D0 and two diffracted harmonics: H0 (the main-line
tail, coupling pi alpha J_0, detuning eta) and Hn (the satellite, coupling
pi alpha J_{-n} e^{∓ i n q t}, detuning eta - n q). Constant coefficients in the
rotating frame -> one 3x3 matrix exponential. The unfilled layer is the exact
2x2 unmodulated Laue matrix over 1 - t. No fit parameter.

Sign convention checked by trying all four (J_{±n}, e^{±inqt}) variants:
(J_{-n}, e^{+inqt}) is best in every case (e.g. S 0.5, alpha 2, n 1:
0.085 vs 0.33–0.71 for the others) — consistent with the kinematic expansion in
`Drive.eta_channel`.

**Result 1 — line shape.** max |I_model - I_TT| / peak over the window
|eta - n q| < q/2, t in {0.31, 0.50, 0.70, 0.90}:

| S, alpha, n | three-mode | satellite only (no H0) |
|---|---|---|
| 1.5, 2, +1 | 0.025–0.105 | 0.062–0.191 |
| 0.5, 2, +1 | 0.053–0.085 | 0.238–0.267 |
| 1.5, 1, -1 | 0.019–0.048 | 0.011–0.056 |
| 3, 2, 2    | 0.017–0.060 | 0.011–0.049 |

The tail mode matters exactly when the satellite is weak relative to the main
line (J_0 large, J_n small: second row, 0.25 -> 0.07). Residual 5–10 % of the
peak remains; candidates are D0's own harmonics and the other satellites'
tails, not included.

**Result 2 — the "peak shift" of Entry 3 is reproduced.** t = 0.3137, S = 1.5,
peak position minus n q (n = +1 / -1):

| alpha_ac | alpha | three-mode | full TT |
|---|---|---|---|
| 20 | 0.5 | +0.201 / +2.857 | +0.222 / +2.828 |
| 20 | 1.0 | +0.135 / +2.869 | +0.137 / +2.849 |
| 40 | 0.5 | +1.548 / +1.145 | +1.444 / +1.234 |
| 40 | 1.0 | +1.474 / +1.154 | +1.354 / +1.252 |
| 80 | 0.5 | -0.553 / -1.061 | -0.562 / -0.925 |
| 80 | 1.0 | -0.539 / -0.986 | -0.553 / -0.837 |

Sign flips with alpha_ac are reproduced; values within ~0.15.

**Reading.**
- The apparent peak shift is now explained quantitatively, not just
  classified: it is interference between the satellite Hn and the main-line
  tail H0 + unfilled-layer tail, whose relative phase e^{i n q t} runs with the
  drive. The same phase produces the acoustic-period micromotion of Entry 2 and
  is the mechanism behind the first-wavelength asymmetry of Entry 8.
- So all four "anomalies" found so far (Entries 2, 3, 8, and the +12 % of
  C-1.10) are one mechanism: coherent overlap of a satellite with the main
  line's off-Bragg tail. They scale as alpha / alpha_ac (tail over
  satellite spacing) and vanish for well-separated channels.
- For the paper: the three-mode model is a compact, closed-form description of
  the whole time-resolved rocking curve near each satellite, good to ~5–10 %
  of the peak. The collapse law (claim 1) is its satellite-only limit.

---

## 2026-09-30 — Entry 15: a ramp through the Bessel ladder — overshoot and a dark window at weak coupling

**Hypothesis.** With a finite rise, the depth profile of the modulation runs
from 0 to S. The local satellite coupling J_n(S A) is not monotonic in S A:
J_1 peaks at 1.841 and changes sign at 3.832. The accumulated-coupling law
(Entry 5) then predicts, even at weak coupling where sin^2 is monotonic in
Theta:
- S beyond the maximum of J_n: while the ramp is inside the crystal, part of
  the depth has *stronger* coupling than the final one -> overshoot.
- S beyond the zero of J_n: positive and negative coupling regions cancel ->
  the satellite appears, goes dark, and reappears.

**Setup.** `explore/nb/e15_cancel.py` and an inline scan. alpha = 0.5 (weak:
pi alpha Theta < pi/2 throughout, so no Pendellösung), alpha_ac = 20,
n = 1, rise = 1 (and 2), centre intensity.

**Result.** Law check: max |I - sin^2(pi alpha Theta_1)| = 0.006–0.013 in all
four runs. Time series (t = 0 to 2 in steps of 0.2; signed Theta below):

| S | I(t) | signed Theta_1(t) | peak / final |
|---|---|---|---|
| 1.5 | 0.00 0.00 0.00 0.03 0.10 0.23 0.38 0.50 0.57 0.59 0.59 | monotone to +0.56 | 1.00 |
| 2.0 | 0.00 0.00 0.01 0.05 0.15 0.30 0.46 0.57 0.62 0.62 0.62 | monotone to +0.58 | 1.00 |
| 3.0 | 0.00 0.00 0.01 0.08 0.19 0.29 0.38 0.38 0.32 0.27 0.26 | up to +0.42, back to +0.34 | **1.52** |
| 5.0 | 0.00 0.00 0.03 0.10 0.10 0.05 0.01 0.02 0.16 0.23 0.24 | +0.20, through 0, to -0.33 | 1.00, with a **dark window** at t ~ 1.2 |

rise = 2, S = 5: peak 0.351 at t = 1.28, final 0.242 (peak/final = 1.45).

**Reading.**
- A 50 % overshoot with no dynamical diffraction involved: at alpha = 0.5
  the satellite is nearly kinematic. The overshoot comes from the drive's own
  Bessel ladder sampled across the ramp. This is a *second, independent*
  overshoot mechanism, distinct from the Pendellösung overshoot (claim 3,
  first maximum of int J_0 = 1.47x). A measurement must separate them: vary
  the rise (this one moves and changes size; the Pendellösung one does not,
  C-1.7(b)) or vary L at fixed S (this one does not depend on alpha).
- The dark window at S = 5 is a clean, qualitative prediction: a satellite
  that switches off for part of the ramp transit and switches back on. Its
  timing is set by where S A(s) = 3.832 sits in the ramp, i.e. by the
  transducer's rise profile — a possible diagnostic of the rise shape itself.
- Caveat: S = 3–5 means u0 ~ 0.1–0.15 A for Si 220 (Entry 6 scaling), strain
  ~2–3e-4 at 2 GHz; large but not absurd. The model has no acoustic
  nonlinearity or attenuation.

---

## 2026-09-30 — Entry 16: validity of the filled-depth law

`explore/nb/e16_validity.py`. S = 1.5, n = 1, centre intensity, 101 times,
max |I - sin^2(pi alpha_1 t)|:

| alpha / alpha_ac | alpha_ac = 5 | 10 | 20 | 40 |
|---|---|---|---|---|
| 0.013–0.025 | — | — | 0.005 (alpha .5) | 0.003 / 0.006 |
| 0.05 | — | 0.011 | 0.010 | 0.013 |
| 0.10 | 0.023 | 0.026 | 0.035 | 0.052 |
| 0.20 | 0.063 | 0.107 | 0.146 | 0.175 |
| 0.40 | 0.312 | 0.341 | 0.429 | — |
| >= 0.8 | 0.61–0.95 | 0.74 | — | — |

**Rule for the paper (revised in Entry 19 — use that one):** alpha / alpha_ac <= 0.1 (equivalently the static
main-line width 2 pi alpha is <= 10 % of the satellite spacing q) keeps the
centre law within 0.05 absolute. At fixed ratio the error grows slowly with
alpha — more Pendellösung turns accumulate the neglected off-resonant phase —
so for many oscillations in transit (alpha_n >~ 3) use <= 0.05.
The law fails completely at alpha / alpha_ac >~ 0.4, where channels overlap
(C-1.6).

**Entry 12 claims table — additions.**

| # | Claim | Evidence | Status |
|---|---|---|---|
| 10 | S from tau_2 / tau_1 | Entry 13 | **V** (simulated measurement, ~1 %, needs alpha_2 >= 0.5) |
| 12 | Three-mode two-layer closed model reproduces I(eta, t) near each satellite to 5–10 % of peak, including apparent shifts | Entry 14 | V |
| 13 | Ramp through the Bessel ladder: 1.5x overshoot at weak coupling (S = 3); dark window (S = 5) | Entry 15 | V |
| 14 | Validity: alpha / alpha_ac <= 0.05 -> centre error <= 0.03; <= 0.1 -> <= 0.08 | Entries 16, 19 (30 random conditions) | V |

---

## 2026-09-30 — Entry 17: Bragg satellite width — narrowing stops at the satellite Darwin width

**Hypothesis.** Satellite n in Bragg is a non-absorbing Bragg slab of coupling
a = pi alpha_n and thickness t (u = (eta - n q)/2):
|u| < a: R = a^2 sinh^2(s t) / (s^2 + a^2 sinh^2(s t)), s = sqrt(a^2 - u^2);
|u| > a: sinh -> sin with s = sqrt(u^2 - a^2).
Thin: t FWHM -> 5.566 (kinematic). Thick: the total-reflection band
|eta - n q| < 2 pi alpha_n, full width 4 pi alpha_n = 12.6 alpha_n.

**Setup.** `explore/nb/e17_bragg_width.py`. alpha_ac = 20, n_z = 1500,
eta step 0.03, 39 times, raw I, same `peak_fwhm` on data and theory.

**Result.**

| S, alpha | alpha_n | max rel. dev. (alpha_n t > 0.3) | FWHM / alpha_n at t = 1, sim / th |
|---|---|---|---|
| 1.5, 1 | 0.56 | 1.4 % | 18.38 / 18.38 |
| 1.5, 2 | 1.12 | 6.1 % | 15.18 / 15.20 |
| 1.5, 4 | 2.23 | 11.5 % (earliest point) | 13.37 / 13.49 |

Example (1.5, 2): alpha_n t = 0.56 / 0.90 / 1.06 gives FWHM 20.1 / 17.9 / 17.2
against theory 20.5 / 17.9 / 17.2 and the kinematic 5.566/t = 11.1 / 6.9 / 5.8.

**Reading.**
- In Bragg the width leaves the kinematic 1/t law early (already 2x at
  alpha_n t ~ 0.5) and **saturates** at the satellite's own Darwin width,
  approaching 13.4 alpha_n (thick-slab value, a little above the 4 pi alpha_n
  band because FWHM includes the edges).
- So in both geometries the width carries the dynamical timescale:
  Laue — narrowing below the finite-depth limit, minimum 0.71x at
  alpha_n t ~ 0.86 (Entry 4); Bragg — narrowing stops at t ~ 1/alpha_n.
- A width measurement in Bragg is therefore a second, independent way to read
  alpha_n (and with it S, via J_n) from the saturated FWHM / 13.4.
- Earliest times (alpha_n t < 0.1): data FWHM ~4 — the window-limited /
  unfilled-crystal fringe problem of C-1.8, as in Laue. Not usable there.

---

## 2026-09-30 — Entry 18: session 1's overshoot table, explained without a fit

C-1.7(a) reported transient max / steady ratios of the channel integrals at
alpha_ac = 20 and called their alpha dependence "weak". The collapse law
predicts ratio = max_{t<=1} F(2 pi alpha_n t) / F(2 pi alpha_n), with
F(x) = int_0^x J_0. The drive amplitude used in session 1 is not written in
C-1.7(a); S = 1 is the value used in C-1.7(b) and reproduces the table:

| alpha | n=1 measured (C-1.7a) | n=1 predicted, S=1 | n=2 measured | n=2 predicted |
|---|---|---|---|---|
| 1 | 1.06 | 1.02 | 1.09 | 1.00 |
| 2 | 2.28 | 2.20 | 1.09 | 1.00 |
| 4 | 1.89 | 1.77 | 1.13 | 1.04 |
| 6 | 1.63 | 1.49 | 1.59 | 1.64 |

(S = 0.5 gives 1.00 / 1.07 / 2.04 / 1.19 for n = 1 and does not fit.)

**Reading.** Every entry agrees to within ~0.1, the residual being the same
sign as the static-tail excess of Entry 2. The "weak, irregular alpha
dependence" was the oscillating function F(2 pi alpha_n) / F_max sampled at a
few alpha: large when the transit ends near a minimum of F (alpha_1 = 0.88 at
alpha = 2, where 2 pi alpha_1 = 5.5 is near the first minimum of F at 5.52),
small otherwise. The n = 2 entries stay ~1 until alpha_2 is large enough to
reach the first maximum inside the transit (alpha = 6). Nothing in that table
is unexplained any more. **Confirmed:** the session-1 script uses
`Drive(S=1.0, alpha_ac=20.0)` (explore/make_feature_figure.py:98).

---

## 2026-09-30 — Entry 19: random conditions — the law is not a special-point artefact

`explore/nb/e19_random.py`, seed 20260930. 30 conditions: S in [0.3, 3.5],
alpha in [0.3, 4], **non-integer** alpha_ac in [20, 45], n in {+1, +2, -1},
rise 0 (half the cases) or U[0, 0.6], 40 random times in [0, 1 + rise].
Law: accumulated coupling (Entry 5). Max abs error of the centre intensity:

| alpha / alpha_ac | cases | worst | median |
|---|---|---|---|
| <= 0.05 | 12 | 0.031 | 0.013 |
| 0.05–0.10 | 13 | 0.075 | 0.048 |
| 0.10–0.20 | 5 | 0.112 | 0.078 |

**Reading.** Holds across random, non-integer, ramped conditions; nothing
depends on integer alpha_ac orthogonality or on special S. The Entry 16 rule
was slightly optimistic. **Revised rule:** alpha / alpha_ac <= 0.05 for
errors <= 0.03; <= 0.1 for <= 0.08.

---

## 2026-10-01 — Entry 20: experiment mapping with verified X0h constants

**Source.** X0h (x-server.gmca.aps.anl.gov), Si, sigma, symmetric Laue,
queried 2026-10-01 via `x0h_form.pl` (parameters in the script). Convention
checked by hand: X0h's Laue extinction length is lambda cos(theta_B) / |chi_rh|
(Si 220, 8 keV: 1.5498 A x 0.91495 / 9.3209e-6 = 15.21 um vs 15.203 reported),
which is the intensity Pendellösung period — the code's Lambda_P
(I = sin^2(pi L / Lambda_P)). Entry 6's guess of 15 um is confirmed.
Shear speed v_T = sqrt(c44 / rho) = sqrt(79.6 GPa / 2329 kg m^-3) = 5846 m/s
(c44 textbook value, not from X0h).

**Key identity.** alpha / alpha_ac = Lambda_ac / Lambda_P. The validity rule
(Entry 19) is a bound on the **acoustic wavelength relative to the
Pendellösung length**, and the minimum drive frequency
f_min = v_T / (r Lambda_P), r = 0.05 or 0.1, does not depend on L.

**Table.** `explore/nb/e20_units_x0h.py`. S = 1.5, L = 3 Lambda_P (alpha = 3,
alpha_1 = 1.67 oscillations in transit). Strain at r = 0.1 is K u0 = S d / Lambda_ac.
Contrast = 1 / cosh(2 pi alpha'' J_1(S)) at t = 1 with alpha''/alpha' = chi_ih/chi_rh.
Transmission = exp(-mu0 L / cos theta_B) (the common factor the code gauges out).

| E keV | hkl | Lambda_P um | L um | t_tr ns | tau_1 ns | f_min GHz (r .05 / .1) | strain (r .1) | kappa | contrast | transmission |
|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 220 | 15.20 | 45.6 | 7.8 | 4.7 | 7.7 / 3.9 | 1.9e-4 | 0.962 | 0.93 | 0.49 |
| 8 | 440 | 15.42 | 46.2 | 7.9 | 4.7 | 7.6 / 3.8 | 9.3e-5 | 0.866 | 0.87 | 0.32 |
| 17.48 | 220 | 36.29 | 109 | 18.6 | 11.1 | 3.2 / 1.6 | 7.9e-5 | 0.965 | 1.00 | 0.85 |
| 17.48 | 440 | 54.34 | 163 | 27.9 | 16.7 | 2.2 / 1.1 | 2.7e-5 | 0.867 | 0.99 | 0.78 |
| 30 | 220 | 63.34 | 190 | 32.5 | 19.4 | 1.9 / 0.9 | 4.5e-5 | 0.966 | 1.00 | 0.95 |
| 30 | 440 | 98.76 | 296 | 50.7 | 30.3 | 1.2 / 0.6 | 1.5e-5 | 0.872 | 1.00 | 0.92 |

(400 rows in the script output; they sit between 220 and 440.)

**Reading.**
- **Correction to Entry 6.** Its 2 GHz plan for Si 220 at 8 keV used
  alpha_ac = max(20, 5 alpha), i.e. r = 0.2 — outside the rule (error up to
  ~0.17). At 8 keV the resolved regime needs 4–8 GHz.
- **Higher energy and higher order are better on every axis.** Lambda_P grows,
  so the required frequency and strain fall, tau_n rises to tens of ns (easy
  timing at a storage ring, even with ~100 ps bunches), and absorption
  disappears (alpha'' / alpha' < 0.4 %). **Recommended poster example: Si 440
  at 17–30 keV, L ~ 150–300 um, shear wave ~1–2 GHz, strain ~2–3e-5.**
- tau_1 = Lambda_P / (v_T J_1(S)) is confirmed thickness-independent in the
  table (L changes, tau_1 tracks Lambda_P only).
- Still unverified: c44 (textbook), whether a clean ~1 GHz shear wave can be
  launched through 0.2–0.3 mm of Si, acoustic attenuation at GHz.

---

## 2026-10-01 — Entry 21: prior art, second pass

**Read in full (2 pages):** DeCamp, Reis, Bucksbaum et al., "Picosecond
Switching of X-rays by Using the Borrmann Effect", APS Activity Report 2001
(companion to Nature 413, 825 (2001)). Facts from the text:
- Ge (001), 280 um, asymmetric 202 Laue, 10 keV; crystal thick against
  absorption, so only the Borrmann mode transmits.
- The disturbance is a **single short strain pulse**, "spatially short compared
  to ... Pendellösung depth", launched by a 50 fs laser.
- Mechanism stated: to switch efficiently "the short acoustic pulse must travel
  into the crystal about one quarter of a Pendellösung depth, or ~400 ps";
  observed switching was faster (supersonic, ambipolar plasma diffusion).
- **No satellites, no periodic drive, no Bessel functions.**

So the main-line "temporal Pendellösung" prior art is a pulse acting as a
moving interface between two perfect crystals. Our n = 0 law (Entry 1) is the
CW-grating version of that picture; the satellite results are outside it.

**Found by search (abstracts only):**
- Polikarpov, phys. stat. sol. (b) 143 (1987) and 149 (1988): dynamical
  diffraction (Laue, Bragg) in crystals weakly deformed by ultrasonic waves —
  the steady-state theory most likely to contain the J_n(S)-scaled satellite
  extinction length. **Must be read before claiming anything about the
  steady-state satellite Pendellösung.**
- Iolin & Entin: suppression of anomalous transmission (Borrmann) by a
  deformation superlattice — steady state, absorption-centred.
- Zolotoyabko, Shilo et al., Rev. Sci. Instrum. 73, 1643 (2002) and related:
  stroboscopic diffraction from 0.58 GHz standing SAWs in LiNbO3 at the APS;
  satellites modulated at 2f — periodic steady state, phase-locked, not a
  formation transient.
- Yoshimura & Hirano, J. Synchrotron Rad. 16, 601 (2009): temporal oscillation
  of Pendellösung fringes without acoustic drive — unrelated.

**Novelty position after this pass.** Nothing found on (i) the transient
formation of satellites behind a travelling front, (ii) the time-domain
collapse onto the static thickness curves, (iii) Bessel-ordered formation,
(iv) the accumulated-coupling law for a ramp, or (v) the period-ratio
calibration of S. The search was not exhaustive (web search, no database
access). The steady-state ingredient (satellite extinction length ~ Lambda/|J_n|)
should be attributed to the ultrasonic-superlattice literature, pending
Polikarpov.

Sources: aps.anl.gov/.../apsar2001/DECAMP1.PDF ;
onlinelibrary.wiley.com/doi/abs/10.1002/pssb.2221490203 ;
onlinelibrary.wiley.com/doi/pdf/10.1002/pssb.2221430102 ;
link.springer.com/article/10.1134/1.1376457 ;
pubs.aip.org/aip/rsi/article-abstract/73/3/1643/1072347 ;
researchgate.net/publication/26773530.

---

## 2026-10-01 — Entry 22: the time-resolved rocking curve as a picture

`explore/make_map_figure.py` -> `figures/poster_map.png`. Laue, n = 1,
S = 1.5, alpha = 3, alpha_ac = 60 (ratio 0.05), 321 eta x 99 t.
(a) full TT, (b) three-mode closed model (Entry 14), (c) centre cut with
sin^2(pi alpha_1 t). max |model - TT| / max TT = 0.065.

What the map shows in one look: the satellite appears broad, narrows as
1/t, blinks off at alpha_1 t = 1 (t = 0.6) with the line splitting into two
side lobes, and returns narrower. The fringes fanning out to both sides are
the finite-depth side lobes. This is the abstract's "time-resolved rocking
curve I(eta, t)" with its analytic counterpart next to it.

**Entry 12 table — status changes.**
- Claim 11 (experimental accessibility): E -> **V for the X-ray constants**
  (X0h, Entry 20); c44 and transducer feasibility remain estimates. The example
  changes to Si 440 at 17–30 keV, ~1–2 GHz (Entry 20).
- Claim 12 (three-mode model): add the map, 6.5 % max deviation at ratio 0.05.

---

## 2026-10-01 — Entry 23: at the dark point the satellite is not gone — fringes flow inward in eta

**Question (from the user's review).** At alpha_n t = 1 the centre is dark.
Does the satellite vanish, or is its intensity redistributed in reciprocal
space? And how do the side lobes move?

**Setup.** `explore/nb/e23_dark.py`. Laue, n = 1, (S, alpha, alpha_ac) =
(1.5, 3, 60) and (1.5, 2, 40) (ratio 0.05). eta step 0.1 over ±2.5 q,
39 times. Tracked: centre I(nq), window integral p (|eta - nq| < q/2), total
over all eta, and the highest local maximum on each side of nq. Theory: the
Laue slab line shape of the satellite, same peak finder.

**Result (case 1.5, 3, 60; alpha_1 = 1.674).**

| alpha_n t | centre (th) | p / (pi^2 alpha_n) (th) | side lobes (eta - nq) t, sim (th) |
|---|---|---|---|
| 0.341 | 0.772 (0.770) | 1.439 (1.452) | — (single peak) |
| 0.425 | 0.946 (0.945) | **1.450** (1.453) max | -8.41 / +8.18 (±8.58) |
| 0.508 | **0.999** (0.999) max | 1.335 (1.331) | -8.16 / +7.91 (±8.40) |
| 0.843 | 0.220 (0.224) | 0.674 (0.677) | -7.09 / +7.02 (±7.26) |
| 0.927 | 0.050 (0.052) | **0.669** (0.684) min | -6.96 / +6.66 (±6.85) |
| 1.010 | **0.002** (0.001) dark | 0.755 (0.774) | -6.74 / +6.42 (±6.36) |
| 1.262 | 0.542 (0.536) | 1.188 (1.198) | -4.20 / +4.17 (±4.23) |
| 1.345 | 0.787 (0.782) | 1.252 (1.263) | -2.95 / +2.84 (±3.06) |
| 1.513 | 0.998 (0.998) | 1.153 (1.179) | -12.35 / +12.05 (±12.18) — next pair |

Second case (1.5, 2, 40) agrees to the same level; at its dark point
(alpha_n t = 1.008): centre 0.001, p = 0.751 pi^2 alpha_n (th 0.771),
lobes -6.46 / +6.28 (th ±6.37).

**Reading.**
1. **The satellite does not vanish at the dark point.** The centre is at
   0.1–0.2 % of its maximum while the window integral is 52 % of its maximum
   (0.755 / 1.450). The intensity sits in two side lobes at
   (eta - nq) t ~ ±6.4. *Formation defined as "growth of the peak intensity" is
   the wrong definition*: the peak blinks, the satellite does not.
2. **The line shape evolves by inward flow of Pendellösung fringes.** The
   maxima of the Laue slab sit where t sqrt(a^2 + u^2) = (k + 1/2) pi, i.e.
   approximately (eta - nq) t = 2 pi sqrt((k + 1/2)^2 - (alpha_n t)^2). As alpha_n t
   grows, each fringe k moves toward nq and merges into the centre at
   alpha_n t = k + 1/2 — **each centre maximum is the arrival of a fringe**,
   and the dark point is the moment between arrivals, when the next fringe
   pair is still at ±6.4 / t. The fringe pair at 1.345 (±3.0) merges at ~1.43,
   and a new pair appears from outside (±12.2 at 1.513). So it is not a
   bifurcation of the central peak: it is a conveyor of fringes.
   The approximate formula gives 6.97 at the dark point vs 6.4–6.7 measured
   (the 1/(a^2 + u^2) envelope shifts the true maxima inward); the exact line
   shape (theory column) matches the simulation to ~5 % in (eta - nq) t.
3. **Integrated and centre clocks are out of phase.** Integral maxima at
   2 pi alpha_n t = zeros of J_0 (alpha_n t = 0.383, 1.377), centre maxima at
   alpha_n t = 0.5, 1.5: the integral leads by ~1/8 of a period, and its
   minimum (0.88–0.93) comes *before* the centre's zero (1.0). The integral
   oscillation is damped (J_0 envelope ~ (alpha_n t)^(-1/2)), the centre's is not.
   Same period, different phase and damping — so "the formation time" depends
   on the observable, and the two must not be mixed.
4. The total over all eta is dominated by the n = 0 channel and the main-line
   tail; it is reported in the script output but carries no satellite-specific
   message.

**Measurement artefact.** For alpha_n t < 0.4 the side-lobe finder returns
values near 0 — it catches the shoulder of the central peak (tail asymmetry)
rather than the first sinc side lobe at ±8.99 / t. Rows before the first
fringe arrival are therefore marked "single peak".

---

## 2026-10-01 — Entry 24: one front, three satellites at different stages

`explore/make_ladder_map.py` -> `figures/poster_ladder_map.png`. Laue, S = 3
(J_1 = 0.339, J_2 = 0.486, J_3 = 0.309), alpha = 3, alpha_ac = 60 (ratio 0.05),
alpha_n = 1.02 / 1.46 / 0.93. I(eta, t) in ±40 around each n q on a shared
colour scale and t axis; centre cuts against sin^2(pi alpha |J_n| t).

At the two marked moments (centre value, window integral over ±40):

| t | n = 1 | n = 2 | n = 3 |
|---|---|---|---|
| 0.344 (n = 2 at its maximum) | 0.788, 13.7 | **0.996**, 17.0 | 0.687, 11.8 |
| 0.684 (n = 2 dark) | 0.674, 8.0 (falling) | **0.000**, 8.8 | 0.807, 8.4 (rising to max) |

**Reading.**
- Under one acoustic front the three orders sit at different stages of the
  same evolution: at t = 0.684 the n = 2 line is split around a dark centre
  while n = 3 is a single narrow line approaching its maximum and n = 1 is
  past its maximum. This is tau_n = Lambda_P / (v_s |J_n(S)|) shown in one
  picture, and it is a stronger statement than the formation-order plot:
  the orders are not just ordered in time, they are in *different line-shape
  phases* simultaneously.
- The dark n = 2 line still holds 52 % of its maximal window integral
  (8.8 / 17.0) — Entry 23 again, now with neighbours as an internal reference.
- n = 3 centre runs ~3 % below its law near its maximum (0.96 vs 1.00);
  n = 3 is closest to the edge of the eta grid used for the other orders'
  tails and is the weakest; not investigated.
- Poster candidate: this figure, or panel d of it, could replace panel c of
  `poster_mapping.png`.

---

## 2026-10-01 — Entry 25: at the Bessel zero the slow clock stops and only the acoustic ripple is left

**Question (user's review).** As S -> j_{1,1} = 3.8317, J_1 -> 0 and
tau_1 = Lambda_P / (v_s |J_1|) -> infinity. Where and how does the
universal law break, and what is left at the zero?

**Setup.** `explore/nb/e25_zero.py` plus an inline fine-time run. Laue, n = 1,
centre intensity, 51 times (fine run: 801).

**Result (a) — S scan, alpha = 3, alpha_ac = 60.**

| S - Z | J_1 | max_t I sim (law) |
|---|---|---|
| -0.30 | +0.124 | 0.845 (0.850) |
| -0.10 | +0.041 | 0.139 (0.141) |
| -0.03 | +0.012 | 1.78e-2 (1.31e-2) |
| -0.01 | +0.004 | 6.95e-3 (1.45e-3) |
| 0 | 0 | 5.59e-3 (0) |
| +0.01 | -0.004 | 6.67e-3 (1.45e-3) |
| +0.03 | -0.012 | 2.01e-2 (1.29e-2) |
| +0.10 | -0.040 | 0.147 (0.135) |

**(b) — residual at the zero, max_t I:**

| alpha | alpha_ac = 30 | 60 | 120 |
|---|---|---|---|
| 1 | 3.2e-3 | 6.3e-4 | 1.8e-4 |
| 2 | 1.3e-2 | 2.5e-3 | 7.1e-4 |
| 3 | 2.9e-2 | 5.6e-3 | 1.6e-3 |

Scales as alpha^2 (1 : 4 : 8.9) and roughly as alpha_ac^-2 (x5.1, x3.5 per
doubling).

**(c)** The minimum of max_t I over S sits at the zero itself (S - Z = 0.0000
for alpha_ac = 30 and 60, step 0.005): no shift of the zero.

**Fine-time run at the zero** (alpha = 3, 801 times): the residual is a ripple
at the acoustic frequency (FFT peak at 59.9 per unit t for alpha_ac = 60;
60.9 for 60.5) with a **constant envelope** (7.1e-3) throughout the transit.
A first reading that "the residual vanishes at t = 1" (I(1) = 4e-6) was a
sampling accident: t = 1 is a whole number of acoustic periods for integer
alpha_ac and lands on a node; with alpha_ac = 60.5, I(1) = 1.2e-3.

**Reading.**
- Away from the zero the law holds to the usual accuracy. It is overtaken by
  the residual when (pi alpha J_1)^2 ~ residual, i.e. |J_1| ~ 0.008 here, or
  **|J_n(S)| alpha_ac <~ 0.5** in general (equivalently J_n L / Lambda_ac <~ 0.5).
- What is left at the zero is not a slow satellite at all: it is the
  acoustic-period micromotion of the main-line tail (Entries 2, 14), with
  amplitude ~ (alpha / alpha_ac)^2 and no growth. **The slow clock stops; the
  fast clock (acoustic period) keeps ticking.** The two timescales separate
  cleanly at the zero.
- Physically interesting by-product: close to the zero, the same crystal
  (alpha = 3, strongly dynamical for the main line) hosts a *kinematic*
  satellite (alpha_1 << 1). S tunes the satellite from dynamical to kinematic
  without changing the crystal — the S-scan rows above go from Pendellösung
  (alpha_1 = 0.37 at S - Z = -0.3, not yet one oscillation) to pure t^2 growth.
- Not to be called critical slowing-down: nothing collective, just
  tau_1 ~ 1 / |S - Z| near a simple zero (J_1 linear in S - Z).

---

## 2026-10-01 — Entry 26: regime map — accuracy and clock are separate knobs

**Idea (user's review).** In the dimensionless model v_s enters only through
t_tr. Two ratios should then control two different things:
x = Lambda_ac / Lambda_P = alpha / alpha_ac (how clean the filled-depth mapping
is) and alpha_n = L |J_n| / Lambda_P (how many satellite Pendellösung turns fit
in the transit, i.e. the clock in units of t_tr).

**Setup.** `explore/make_regime_map.py` -> `figures/regime_map.png`. Laue,
n = 1, S = 1.5 fixed; 7 x 7 log grid x in [0.0125, 0.4], alpha_1 in [0.25, 4];
alpha = alpha_1 / J_1, alpha_ac = alpha / x; 121 times; colour = max_t of
|I(q, t) - sin^2(pi alpha_1 t)|.

**Result.** Error by column (x), all seven alpha_1 rows:

| x | 0.0125 | 0.022 | 0.040 | 0.071 | 0.126 | 0.22 | 0.40 |
|---|---|---|---|---|---|---|---|
| range over alpha_1 | 0.005–0.013 | 0.007–0.019 | 0.012–0.035 | 0.015–0.062 | 0.030–0.113 | 0.050–0.193 | 0.18–0.40 |
| typical | 0.01 | 0.02 | 0.03 | 0.06 | 0.10 | 0.18 | 0.36 |

**Reading.**
- **The error is set by x alone**, to a good approximation
  max error ~ 0.9 x = 0.9 Lambda_ac / Lambda_P, with no systematic trend in
  alpha_1 from 0.25 to 4 turns. The two knobs separate: x fixes how faithful the
  mapping is; alpha_n fixes what the satellite does in time.
- Physically: at fixed material (v_s fixed) the drive frequency moves x only —
  higher f, cleaner mapping — while the clock tau_n = Lambda_P / (v_s |J_n|)
  does not depend on f at all. The thickness L moves alpha_n only. S moves
  alpha_n only (through J_n). So an experiment can tune cleanliness and clock
  independently, which promotes the design rule of Entry 20 to a regime
  statement.
- Scatter below the typical value comes from near-integer alpha_ac cells
  (e.g. alpha_1 = 0.25, x = 0.22 has alpha_ac = 2.0: error 0.050 instead of
  ~0.18) — the finite-slab orthogonality of integer alpha_ac suppresses the
  tail at the satellite centre. Real crystals will not sit on integers; quote
  the typical value.
- This revises Entry 16's remark that the error "grows slowly with alpha at
  fixed ratio": at fixed S it does not; the Entry 16 trend came from
  alpha_1 and the integer-alpha_ac structure varying together.
- Caveat: stroboscopic model. Any genuine dependence on v_s beyond t_tr
  (e.g. the X-ray response time, quasi-static breakdown, handoff 7.5) is outside
  this map by construction.

**Entry 21 addendum (user, 2026-10-01).** Time-resolved rocking curves of
ultrasonically excited crystals have been measured for a long time (e.g. Liss
et al., 1998: liss.freeshell.org/kdl/papers/Liss98-2.pdf), and the time
dependence of SAW satellites is known. Those works look at the phase of an
already established periodic (often standing) acoustic field. The question
here — the line shape of a satellite being *created* behind a travelling
front — is a different problem. Liss 1998 not yet read by me.

---

## 2026-10-01 — Entry 27: claims table, additions from Entries 20–26

| # | Claim | Evidence | Status |
|---|---|---|---|
| 15 | At the dark point the satellite keeps ~52 % of its maximal integrated intensity, in two side lobes: formation is not peak growth | Entry 23, test `test_dark_centre_keeps_integrated_satellite` | V |
| 16 | The line shape evolves by inward flow of Pendellösung fringes in eta; each centre maximum is a fringe arrival; lobe positions collapse as (eta - nq) t = G(alpha_n t) | Entry 23 | V (~5 % in (eta - nq) t) |
| 17 | Integrated and centre intensities share the period but not the phase or damping (integral leads by ~1/8 period, damped as J_0) | Entry 23 | V |
| 18 | Under one front, different orders are in different line-shape phases at the same time (n = 2 dark while n = 3 peaks, S = 3) | Entry 24, `figures/poster_ladder_map.png` | V |
| 19 | At a Bessel zero the slow clock stops; only an acoustic-period ripple ~ (alpha/alpha_ac)^2 remains; the law fails for abs(J_n) alpha_ac <~ 0.5 | Entry 25 | V |
| 20 | Mapping accuracy ~ 0.9 Lambda_ac / Lambda_P, independent of the number of turns; f tunes accuracy, L and S tune the clock | Entry 26, `figures/regime_map.png` | V |
| 11 | Experimental accessibility | Entry 20 | V for X-ray constants (X0h); c44, transducer: E |

**Sharper one-sentence story (replaces Entry 12's).** Behind a travelling
acoustic front each satellite is a crystal whose thickness grows at v_s and
whose extinction length is Lambda_P / abs(J_n(S)); its line shape is therefore a
conveyor of Pendellösung fringes flowing toward the satellite centre, so the
peak blinks while the integrated satellite persists, and different orders
under the same front are caught at different stages of the same evolution.

---

## 2026-10-01 — Entry 28: a transient-only satellite, and the transiently dominant forbidden order

**Framework (new, `explore/nb/memory_drive.py`).** At time t the crystal holds
the drive history s = t - zeta over the last t_tr: it is a spatial record of
the transducer signal. Satellite n reads its Bessel-weighted integral
Theta_n(t) = int_0^1 J_n(S A(t - zeta)) sign_n dzeta (sign_n = (-1)^n where the
carrier is flipped by pi). Readout at exact resonance:
- Laue: sin^2(pi alpha |Theta_n|) — rotations about one axis add.
- **Bragg: tanh^2(pi alpha |Theta_n|), exact** — the Riccati equation at
  resonance, dR/dzeta = -i kappa(zeta) (1 + R^2), is solved by
  R = -i tanh(int kappa) for any real depth-dependent kappa.
So Theta_n is geometry-independent; geometry only picks the readout function.
`HistoryDrive` subclasses `Drive` with an arbitrary envelope and an optional
pi carrier flip (no core change).

**Setup.** `explore/nb/e28_transient_only.py`. S_f = j_{1,1} = 3.8317 (J_1 = 0),
linear ramp of length T in {0.25, 0.5, 1, 2, 4}, alpha_ac = 40, alpha in
{0.5, 2}, Laue and Bragg, centre intensities of n = 1–3 (n = 1–4 in a second run),
161 times to t = 1 + T + 0.3.

**Result — transient-only satellite.** n = 1, Laue (Bragg in brackets):

| alpha | T | max_t I_1 (law) | FWHM lifetime | final I_1 | law max err |
|---|---|---|---|---|---|
| 0.5 | 0.25 | 0.023 (0.021) [0.023] | 0.92 | 5e-10 | 0.002 |
| 0.5 | 1 | 0.302 (0.296) [0.274] | 0.73 | 5e-10 | 0.007 |
| 0.5 | 4 | 0.606 (0.603) [0.506] | 2.12 [2.22] | 5e-10 | 0.003 |
| 2 | 0.5 | 0.855 (0.833) [0.686] | 0.92 | 2e-6 | 0.029 |
| 2 | 4 | 0.998 (1.000) [0.997] | 3.45 [3.35] | 2e-6 | 0.010 |

I_1 = 0 before the front, > 0.5 at its peak, and back to 1e-10–1e-6 once the
ramp has left the crystal (Theta_1 -> J_1(S_f) = 0). Max intensity grows with T
and saturates at the quasi-static value sin^2(pi alpha max J_1) = 0.62
(alpha = 0.5); lifetime grows with T for T >~ 1.

**Result — forbidden order transiently dominant.** Max over time vs final:

| geometry, alpha, T | max_t I_n, n = 1 / 2 / 3 / 4 | final I_n |
|---|---|---|
| Laue, 0.5, 4 | **0.606** / 0.464 / 0.370 / 0.140 | **0** / 0.349 / 0.370 / 0.139 |
| Laue, 0.5, 2 | **0.535** / 0.421 / 0.370 / 0.140 | **0** / 0.349 / 0.370 / 0.139 |
| Bragg, 0.5, 4 | **0.506** / 0.401 / 0.330 / 0.134 | **0** / 0.313 / 0.330 / 0.133 |
| Laue, 2, 1 | 0.999 / 0.997 / 0.978 / 0.916 | 0 / 0.329 / 0.226 / 0.914 |

At weak coupling the order that is absent from the final pattern is the
strongest of all orders during the transient (n = 1 leads in ~60 % of the
samples for T = 4). At alpha = 2 every order saturates near 1 and the
comparison stops discriminating.

**Reading.**
- **The steady-state pattern does not predict the transient's dominant
  channel.** A satellite forbidden in steady state (Bessel zero) can be the
  brightest line while the ramp is inside the crystal.
- The effect is identical in Laue and Bragg (law error <= 0.03 in both): it is
  a property of the time-dependent acoustic coupling Theta_n, not of
  Pendellösung. It also exists at alpha = 0.5, where nothing oscillates.
- Memory parameter: T vs t_tr. For T >> 1 the crystal is nearly uniform and
  I_1 tracks the instantaneous J_1(S A(t)) (adiabatic); for T << 1 the
  satellite is a short flash whose size ~ (pi alpha T <J_1>)^2.
- Caveat: S = 3.83 is a large drive (Entry 20 scaling: u0 ~ 0.6 A for
  Si 440). Acoustic nonlinearity is outside the model.

---

## 2026-10-01 — Entry 29: Polikarpov & Skadorov read in full (supplied by the user)

**[P87]** I. V. Polikarpov, V. V. Skadorov, "X-Ray Dynamical Diffraction in
Crystals Being Subjected to the Action of Ultrasonic Waves", phys. stat. sol.
(b) 143, 11 (1987). **[P88]** same authors, "X-Ray Bragg Dynamical Diffraction
in Crystals Weakly Deformed by Ultrasonic Waves", phys. stat. sol. (b) 149,
435 (1988).

What they contain:
- Quasi-static, multi-wave (Floquet) treatment: the non-stationary two-wave
  problem becomes multi-wave diffraction of waves k + p kappa by a stationary
  lattice + ultrasonic superlattice. The coupling coefficients are written with
  Bessel functions explicitly, chi_{p,p'} = (-1)^{p-p'} e^{i(p'-p)phi}
  J_{p'-p}(tau a) C chi_tau (P87 eq. 1) — our J_n(S), S = h.u0 = tau a.
- They then expand in xi = (1/2) tau a << 1 (weak deformation) and keep only
  p = 0, ±1 (six waves). The satellite extinction length is
  Lambda_a = (k sqrt((Delta lambda)^2 + 4 xi^2 eta_a))^-1 (P87 eq. 13): ~ 1/xi,
  the small-S limit of our Lambda_P / |J_1(S)|.
- Main message of P87: the non-resonant waves must be kept; they make the
  extinction lengths of different branches differ and give "an intricate
  oscillating dependence of the X-ray intensity on the ultrasonic wave
  amplitude" — steady state. P88: Bragg geometry, two- vs four-level
  resonance depending on alpha, steady state.
- Uniform, fully modulated crystal throughout. No travelling front, no
  transient, no time-dependent amplitude, no rise. S is perturbative, so the
  Bessel zeros are outside their scope.

**Consequences for our claims.**
- The steady-state satellite extinction length (weak-S form) is prior art:
  P87, with earlier Soviet work it cites (Entin 1977/1979; Iolin & Entin 1983;
  Dzyublik 1983; Assur & Entin 1982). Cite them for the ingredient.
- Their emphasis on non-resonant waves is the steady-state counterpart of our
  three-mode corrections (Entry 14); cite as such.
- Still not found anywhere: the transient formation behind a travelling front
  and everything built on it (collapse, fringe flow, Bessel-ordered stages,
  accumulated-coupling law, transient-only and transiently dominant orders,
  history dependence). Our use of the full J_n(S) at arbitrary S (zeros
  included) also goes beyond their xi << 1 truncation.

---

## 2026-10-01 — Entry 30: formation remembers the path

**Setup.** `explore/nb/e30_history.py`; figure `explore/make_memory_figure.py`
-> `figures/poster_memory.png`. Laue and Bragg, alpha = 0.5 (no Pendellösung),
alpha_ac = 40, centre intensities. Law: Theta_n read by sin^2 / tanh^2
(Entry 28). Law error <= 0.008 in every run below, both geometries.

**(a) Same transducer amplitude, different pattern (triangular pulse,
S_max = 3, rise = fall = T).** I_1 on the way up / down at equal S A(t):

| geometry | T | A = 0.3 | A = 0.6 | A = 0.9 |
|---|---|---|---|---|
| Laue | 0.25 | 0.001 / 0.091 | 0.008 / 0.060 | 0.024 / 0.034 |
| Laue | 1 | **0.010 / 0.494** | 0.117 / 0.518 | 0.322 / 0.429 |
| Laue | 4 | 0.146 / 0.510 | 0.560 / 0.578 | 0.544 / 0.355 |
| Bragg | 1 | 0.010 / 0.428 | 0.111 / 0.442 | 0.286 / 0.372 |

n = 2 behaves the same way (Laue, T = 1, A = 0.3: 0.000 / 0.289).
Reading: the satellite reads the drive averaged over the last t_tr, i.e. the
transducer signal delayed by ~t_tr / 2 and smoothed; plotted against the
instantaneous amplitude the response is a loop (figure panel b). Even a slow
pulse (T = 4 t_tr) keeps a large loop: what matters is dA/dt x t_tr, not
T >> t_tr alone. At A = 0.9 for T = 4 the order reverses (up > down) because the
path crosses the J_1 maximum (S A = 1.84) inside the crystal. **Not**
thermodynamic hysteresis — a finite, exactly known memory of length t_tr.

**(b) Satellite after the drive has stopped.** Transducer amplitude 0
(t = 2T + ...): Laue I_1 = 0.105 (T = 0.25), 0.063 (T = 1): the packet is still
in the crystal. It vanishes once the packet has left (T = 4: 0.005, sampled late).

**(c) Multiple revival, slow triangle to S = 10, T = 3, n = 1.** Bright maxima
0.369, 0.105, 0.095, 0.107, 0.368 at t = 1.05, 2.1, 3.5, 4.85, 5.9; dark minima
(I < 0.02, essentially 0) at t = 1.7, 2.63, 4.38, 5.31 — two J_1 zeros crossed
on the way up and again on the way down. Four dark windows at alpha = 0.5, where
the step-onset satellite would never pass through zero: this is the Bessel
ladder of the drive, not Pendellösung. With a fast pulse (T = 0.5) the crystal
averages over the whole ladder and I_1 stays below 0.05 throughout.

**(d) Phase-coded packet: box of length 0.8 with a pi carrier flip at 0.4, S = 1.5.**

| t | 0.2 | 0.4 | 0.6 | 0.8 | 0.9 | 1.0 | 1.2 | 1.4 | 1.6 | 1.9 |
|---|---|---|---|---|---|---|---|---|---|---|
| n = 1 (Laue) | 0.031 | 0.117 | 0.030 | **0.000** | **0.000** | **0.000** | 0.031 | **0.118** | 0.030 | 0.000 |
| n = 2 (Laue) | 0.005 | 0.022 | 0.048 | **0.083** | 0.083 | 0.084 | 0.048 | 0.020 | 0.005 | 0.000 |

Bragg identical to within 0.006. While both halves are inside, the odd order
cancels exactly (the flipped half contributes (-1)^n J_n) and the even order is
at its maximum; when the leading half leaves, n = 1 **revives** to its earlier
peak. A trailing edge erasing what the leading edge wrote, then the erasure
being undone.

**Test.** `tests/test_formation.py::test_history_law_transient_only_and_phase_coded_packet`
(39 tests pass).

**Story candidate (user's phrasing, now supported):** *formation dynamics
remembers the path, not just the instantaneous strain* — with the memory
length equal to the transit time and the memory kernel known in closed form.

**Limits.** Phase jumps other than 0 / pi are not covered by the law (Entry 28
framework). Real transducers: arbitrary A(s) with a pi phase flip is a
standard RF burst technique, but its fidelity inside the crystal (dispersion,
attenuation) is not modelled. S up to 10 is a very large drive.

---

## 2026-10-01 — Entry 31: the brightest order switches in time, lagging the drive by t_tr / 2 — and a grid-rule bug

**Setup.** `explore/nb/e31_switch.py`. Linear ramp to S = 6 over T (4 and 1),
alpha = 0.5, alpha_ac = 40, Laue, centre intensity of n = 1..5; "brightest order"
= argmax_n (ignored while all < 0.02). Compared with the history law and with
the instantaneous ladder sin^2(pi alpha J_n(S A(t))) (no memory).

**Result (converged grid).** Law error 0.005. Switch times (from -> to):

| T | full TT | history law | instantaneous J_n(S A(t)) |
|---|---|---|---|
| 4 | 1->2 at 2.26, 2->3 at 3.01, 3->4 at 3.74, 4->5 at 4.83 | identical | 1.76, 2.52, 3.27, 4.00 |
| 1 | 0.75–0.76 (with one flicker), 0.97, 1.20, 1.91 | 0.76, 0.97, 1.20, 1.91 | 0.45, 0.63, 0.82, 1.00 |

**Reading.**
- Yes: the brightest satellite changes identity n = 1 -> 2 -> 3 -> 4 -> 5 during
  the ramp, at weak coupling, in the full simulation.
- For the slow ramp the switches lag the instantaneous Bessel ladder by
  0.50, 0.49, 0.47 t_tr: when a linear ramp covers the whole crystal, the depth
  average of S A(t - zeta) is S A(t - 1/2) exactly, so the crystal reports the
  drive **half a transit late**. The last switch (4 -> 5) lags more (0.83)
  because the transducer stops ramping at t = 4 while the crystal is still
  filling. For the fast ramp (T = 1) the lag is 0.31–0.38 and the order
  sequence is compressed.
- A clean experimental observable: the sequence of brightest orders is set by
  the drive's Bessel ladder, the delay of each switch by t_tr / 2.

**Bug found and fixed while doing this.** A first run (n_z = 2400) gave law
error 0.09 at n = 5 and missed the last switch; the error *grew* with alpha_ac,
which ruled out the alpha/alpha_ac correction. Convergence: max error at n = 5
was 0.089 / 0.0031 / 0.0004 / 0.0004 for n_z = 2400 / 4800 / 9600 / 19200.
Cause: `core.required_nz` counted only alpha and alpha_ac, but e^{-i phi}
sweeps up to S q per unit depth and e^{i eta zeta} turns eta / 2 pi times — at
S = 6 read at eta = 5 q that is ~11 alpha_ac cycles, not alpha_ac.
**Fix:** `required_nz(alpha, alpha_ac, per_cycle, S=0, eta_max=0)` now uses
max(alpha, alpha_ac, eta_max / 2 pi + S alpha_ac); `cube.*` and
`sweep.sweep_point` pass S and max |eta|. Test
`test_analytic.py::test_required_nz_resolves_large_S_high_channel`; 40 pass.

**Impact on earlier entries (re-checked).**
- Entry 28 (S = 3.83, n_z = 2400): n = 1 unchanged (max 0.605, final 0); n = 3, 4
  move slightly (max 0.370 -> 0.376, 0.139 -> 0.153). Conclusion unchanged:
  the forbidden n = 1 is still the brightest (0.605 > 0.462).
- Entry 30(c) (S = 10, n = 1): law error 0.0073 at both n_z = 2400 and 9600 —
  converged (low channel).
- Entries with S <= 3 and n <= 3 used n_z >= 40–60 alpha_ac by hand, which meets
  the new rule within a factor ~1.5; not re-run.
- Session 1's sweep (C-1.5) used the old default at S up to 3, windows to ±3 q.
  Its t_f values were already judged uninformative; if they are ever used,
  re-run with the new rule.

---

## 2026-10-01 — Entry 32: two co-propagating waves — combination satellites need overlap

**Setup.** `explore/nb/e32_two_waves.py`. phi = S1 A1(s) sin(q1 s) +
S2 A2(s) sin(q2 s), s = t - zeta, both step envelopes, the second launched
tau = 0.4 later; alpha_ac1 = 40, alpha_ac2 = 57, S1 = S2 = 1.2, alpha = 0.5,
Laue, n_z from the extended rule (4676). Channel (n, m) at eta = n q1 + m q2.
Law: the product of the two Jacobi–Anger expansions gives coupling
J_n(S1 A1) J_m(S2 A2), so Theta_nm = int_0^1 J_n(S1 A1(t - zeta)) J_m(S2 A2(t - zeta)) dzeta.

**Result.** Law max error 0.003–0.006 on (1,0), (0,1), (1,1), (1,-1), (2,-1).
Selected time series (centre intensity):

| t | (1,0) | (0,1) | (1,1) | (1,-1) |
|---|---|---|---|---|
| 0.30 | 0.055 | 0.000 | 0.000 | 0.000 |
| 0.60 | 0.169 | 0.012 | 0.006 | 0.007 |
| 1.05 | **0.333** | 0.113 | 0.062 | 0.063 |
| 1.50 | 0.251 | 0.251 | 0.145 | 0.145 |

Time to reach 20 % of each channel's maximum: (1,0) 0.35; (0,1), (1,1), (1,-1)
0.85; (2,-1) 0.80.

**Reading.**
- The history law extends to two waves without change of form.
- Combination satellites exist only in the depth where both waves are present:
  they start growing at t = tau and their clock is the overlap, not the transit.
- Cross-talk: while the delayed wave fills the crystal, the first wave's own
  satellite (1,0) *dims* (0.333 -> 0.251), because in the overlap its coupling is
  J_1(S1) J_0(S2), J_0(1.2) = 0.67. A second acoustic wave switches the first
  wave's satellite down with a delay set by tau.
- No anomalously large transient combination satellite was found in this
  configuration; the user's "transient combination satellites" idea gives
  ordered, overlap-gated formation rather than an anomaly. Scope note: two-wave
  physics widens the study a lot; parked here.

---

## 2026-10-01 — Entry 33: claims table, additions from Entries 28–32, and the new story

| # | Claim | Evidence | Status |
|---|---|---|---|
| 21 | Theta_n = int_0^1 J_n(S A(t - zeta)) sign_n dzeta is geometry-independent; Laue reads sin^2, Bragg reads tanh^2 (exact at resonance for any depth profile) | Entries 28, 30 (law error <= 0.03 / 0.008) | V |
| 22 | Transient-only satellite: forbidden before and after (J_1(S_f) = 0), bright while a ramp is inside | Entry 28, test | V |
| 23 | Steady-state-forbidden order is the brightest of all during the transient (weak coupling) | Entry 28 | V |
| 24 | History dependence: same transducer amplitude, up vs down a pulse, up to 50x different (I_1 0.010 vs 0.494) | Entry 30a, `figures/poster_memory.png` b | V |
| 25 | Satellite persists after the transducer is silent, while the packet is inside | Entry 30b | V |
| 26 | Multiple bright/dark revivals at weak coupling from the drive's Bessel ladder | Entry 30c | V |
| 27 | Phase-coded packet: pi flip cancels odd orders while both halves are inside, then they revive; even orders add | Entry 30d, test | V |
| 28 | The brightest order switches 1 -> 2 -> 3 -> 4 -> 5 during a ramp, lagging the instantaneous ladder by t_tr / 2 | Entry 31 | V |
| 29 | Two waves: combination satellites are overlap-gated; a second wave dims the first wave's satellite by J_0(S2) | Entry 32 | V |
| — | Steady-state satellite extinction length ~ 1/J_1 (weak S) | Polikarpov & Skadorov 1987/1988 and refs (Entry 29) | prior art — cite |

**Story, version 3.** A travelling acoustic front turns the crystal into a
record of the last t_tr of the transducer signal. Each satellite reads that
record through a Bessel kernel, Theta_n = int J_n(S A) dzeta, and the geometry
only chooses how the reading is displayed — sin^2 (Laue) or tanh^2 (Bragg).
Hence satellite formation remembers the path, not just the instantaneous
strain: satellites can exist only in the transient, steady-forbidden orders
can dominate it, equal strains give different patterns on the way up and down,
and a pi-coded packet can erase and restore a satellite. The step-onset
results (temporal Pendellösung, collapse, fringe flow) are the special case
of a record that is a single step.

---

## 2026-10-01 — Entry 34: is the transient-only satellite real? Robustness and honest framing

User's question: should the transient-only peak be claimed as existing?
Quick checks (inline run, Laue, alpha = 0.5, alpha_ac = 40, ramp T = 2,
n_z from the extended rule), final / transient maximum of I_1:

| perturbation | transient max | final | final / max |
|---|---|---|---|
| none (S_f = j_11) | 0.537 | 0.0000 | 0.000 |
| S_f off by ±2 % | 0.541 / 0.534 | 0.0024 / 0.0023 | 0.004 |
| S_f off by ±5 % | 0.546 / 0.528 | 0.0152 / 0.0138 | 0.03 |
| acoustic amplitude loss 5 % across L | 0.533 | 0.0039 | 0.007 |
| loss 10 % | 0.528 | 0.0161 | 0.03 |
| loss 20 % | 0.517 | 0.0663 | 0.13 |

(attenuation modelled as S e^{-gamma zeta} with S at the zero at the entrance.)

**Verdict.**
- Within the model it is real and robust: full TT, law-exact, both
  geometries, and a > 30:1 transient-to-final contrast survives a 5 % error in S
  or 10 % acoustic loss across the crystal. S can be set to ~1 % with the
  period-ratio calibration (Entry 13).
- Its novelty must be framed carefully. For a slow ramp (T >> t_tr) it is just
  the steady-state Bessel ladder sampled on the way to a zero — any amplitude
  sweep shows n = 1 at intermediate S. What the transient adds is (i) the
  t_tr memory (the satellite reports the drive t_tr / 2 late and smoothed),
  (ii) the fast-ramp regime (T <~ t_tr) where it is a flash set by the ramp
  profile inside the crystal, and (iii) that a step onset gives *exactly zero*
  at all times — so the peak is a signature of the amplitude gradient inside the
  crystal, not of the final state.
- Claim it as: "the steady-state pattern does not determine the transient
  one: an order forbidden at the final amplitude can be the brightest line
  during formation", with the slow-ramp limit stated as the trivial
  (adiabatic) end of it. Do not claim it as a new kind of peak.
- Unmodelled risks for an experiment: acoustic nonlinearity at S ~ 4 (strain
  ~ S d / Lambda_ac ~ 7e-5 for Si 440 at ~1 GHz — modest), depth-inhomogeneous
  amplitude beyond simple attenuation, transducer ringing.

---

## 2026-10-01 — Entry 35: robustness and failure boundary of the transient-only satellite

**Hypothesis.** The transient-only character (n = 1 bright during formation,
dark at the end, S_f = j_11) survives realistic departures from the ideal
envelope A(t - zeta) as long as the *final accumulated coupling*
Theta_1(final) stays small compared with its transient maximum. Different
non-idealities should fall into two classes: those that raise the final
background (break Theta_1(final) = 0) and those that only reshape the transient.

**Setup.** New module `explore/nb/nonideal.py` (`GeneralDrive(Drive)`; no core
change): phase phi = S A(s, zeta) sin(theta) + m2(s, zeta) sin(2 theta + psi2),
theta = q (zeta - t), s = t - zeta. Models (one at a time):
1. attenuation A -> A exp(-gamma zeta); reported as amplitude loss 1 - exp(-gamma) across L;
2. dispersion: envelope convolved with a Gaussian of width sig1 sqrt(zeta)
   (diffusive spreading; the precursor ahead of the front is kept — a non-causal
   tail, acknowledged);
3. nonlinearity: second harmonic m2 = eps2 zeta S A^2 (grows with distance and
   with amplitude squared — minimal steepening model), psi2 = 0 and pi/2;
4. ringing after the ramp: R exp(-x/tau) sin(2 pi f x + psi), x = s - T;
5. combinations.
Law: Theta_n(t) = int_0^1 c_{-n}(S A, m2) dzeta, with c the Fourier coefficient
of exp(-i(a sin + b sin 2)); it reduces to J_n when m2 = 0. Script
`explore/nb/e35_robust.py` (sections base / atten / disp / nl / ring / combo /
retune / integ); figure `explore/make_robust_figure.py` ->
`figures/robust_contrast.png`. Laue, n = 1, alpha_ac = 40, n_z from the extended
rule (3865), 81 transient samples + 16 over the final acoustic period
(I_fin = their mean, which avoids the ripple nodes of Entry 25).
Cases: FAST (alpha 1, T 0.5 — the genuine transient, T < t_tr), MID (alpha 0.5,
T 2 — Entry 28), STEP (alpha 1, T 0 — control), SLOW (alpha 0.5, T 8 — adiabatic
reference, not a transient effect).

**Baseline and controls.**

| case | I_max | I_fin | contrast | law err |
|---|---|---|---|---|
| FAST | 0.290 | 3.1e-8 | 9.4e6 | 0.014 |
| MID | 0.539 | 4.8e-10 | 1.1e9 | 0.006 |
| STEP (control) | **0.0018** | 3.1e-8 | — | 0.002 |
| SLOW | 0.621 | 4.8e-10 | 1.3e9 | 0.002 |
| reference S_f = 2.5 (J_1 = 0.50), FAST | 1.000 | 1.000 | 1.0 | 0.008 |

Convergence: doubling n_z leaves I_max and I_fin unchanged to 5 digits (FAST,
STEP). The step control is not exactly zero during the transit: 0.0018 is the
acoustic-period tail ghost ~ (alpha/alpha_ac)^2 of Entry 25 — the floor any
"transient peak" must exceed.

**1. Attenuation (contrast).**

| loss | 1 % | 3 % | 5 % | 10 % | 20 % | 40 % | 60 % |
|---|---|---|---|---|---|---|---|
| FAST | 511 | 59 | 22 | 6.1 | 2.0 | 1.0 | 1.0 |
| MID | 3640 | 395 | 139 | 33 | 7.9 | 2.1 | 1.1 |
| SLOW | 4200 | 457 | 162 | 39 | 9.4 | 2.5 | 1.4 |

I_max barely changes (FAST even grows, 0.29 -> 0.49 at 20 %); only I_fin rises.
Law error 0.014–0.016 throughout. STEP control: I_max ~ I_fin at every loss —
no transient-only peak without a ramp.

**2. Dispersion.** I_fin stays at the floor (3e-8) for every sig1 up to 0.8:
behind the packet A -> 1 and J_1(S_f) = 0 is untouched. The transient peak grows
(FAST 0.29 -> 0.99 at sig1 = 0.8). **The STEP control becomes bright too**
(0.0018 -> 0.035 / 0.12 / 0.34 / 0.97 at sig1 = 0.05 / 0.1 / 0.2 / 0.8): dispersion
turns the step into a ramp inside the crystal. Law error <= 0.016.

**3. Second harmonic (contrast).**

| eps2 | 0.005 | 0.01 | 0.02 | 0.05 | 0.1 | 0.2 |
|---|---|---|---|---|---|---|
| FAST | 7790 | 1920 | 490 | 86 | 25 | 8.7 |
| MID | 54900 | 13700 | 3460 | 577 | 156 | 47 |

I_fin ~ eps2^2 (3.8e-5 / 1.5e-4 / 6.2e-4 for 0.005 / 0.01 / 0.02). psi2 = pi/2
gives the same background within 10 % (TT only). With the correct harmonic sign
(c_{-1}, see corrections) the law holds to 0.015 up to eps2 = 0.2.

**4. Ringing.** In a fixed measurement window (t = 1 + T + 0.5) the background
rises only when the ringing outlives the transit (tau >= 1: contrast ~2000–2300
FAST, ~15000–17000 MID); short (tau = 0.05–0.2) or fast (f = 10) ringing is
averaged away by the depth integral (contrast 1e5–1e6). Once the ringing has died
(t = 1 + T + 0.5 + 6 tau) the background returns to the floor in every case.
I_max is unaffected (0.289–0.315).

**5. Combined.** FAST: atten 10 % + ringing R 0.05 -> 5.9; dispersion 0.1 + NL
0.02 -> 581; atten 5 % + NL 0.02 + ringing 0.02 -> 15. MID: 33 / 3460 / 96.
Combinations behave like the sum of their Theta_1(final) contributions; the law
holds to 0.015.

**6. Essence test — local Bessel zero or depth integral?** With attenuation,
retune S so that int_0^1 J_1(S exp(-gamma zeta)) dzeta = 0 (S* = 4.039 / 4.284 /
4.958 for 10 / 20 / 40 % loss). Contrast, FAST:

| loss | S = j_11 | S* (fraction of depth with abs(J_1) > 0.05) |
|---|---|---|
| 10 % | 6.1 | **1.1e5** (38 %) |
| 20 % | 2.0 | **2.4e4** (71 %) |
| 40 % | 1.0 | **3.8e3** (87 %) |

The satellite is dark at the end although J_1 is far from zero over most of the
crystal: positive and negative coupling regions cancel in the depth integral.
STEP with S* is also dark at the end (I_max 0.0065–0.089, from the moving front).

**7. Integrated satellite (window abs(eta - q) < q/2), FAST.** Ideal: integrated
contrast 119 (centre 9.4e6); atten 20 %: 2.5 (centre 2.0); NL 0.05: 67 (86);
ringing R 0.1: 120 (1e5). The integrated contrast tracks the centre one when the
satellite background dominates, and is capped at ~120 otherwise by the
main-line tail inside the window (Entry 2) — not by the satellite. The
transient-only character is not a centre-intensity artefact.

**Reading.**
- **One control parameter.** Every non-ideality collapses onto
  r = abs(Theta_1(final)) / max_t abs(Theta_1(t)), with contrast ~ 1/r^2 over six
  decades (figure panel a). **30:1 needs r <= 0.18; 10:1 needs r <= 0.32.**
  Example: FAST, 3 % loss: r = 0.13, predicted 59, measured 58.
- **Two classes.** Attenuation, second harmonic and long ringing *raise the
  final background* (they break Theta_1(final) = 0); none of them destroys the
  transient peak (it stays or grows). Dispersion leaves the background alone and
  *reshapes* the transient — it even creates one from a step. Nothing tried
  destroys the transient peak itself.
- **Failure boundaries in physical terms** (Theta_1(final) ~ 0.78 x loss and
  ~ 0.38 eps2 for small values): 30:1 holds for loss <= 0.23 max Theta_1 and
  eps2 <= 0.47 max Theta_1. FAST (max Theta_1 = 0.18): loss <~ 4 %, eps2 <~ 0.09.
  MID (0.52): loss <~ 12 %, eps2 <~ 0.25. A larger transient Theta (longer ramp,
  larger alpha) buys robustness — but the longer the ramp, the closer to the
  trivial adiabatic scan.
- **The essence is the accumulated coupling, not the Bessel zero.** The Bessel
  zero is one way to make Theta_1(final) = 0; spatial cancellation is another
  (section 6). A transient Theta_1 != 0 needs a spatially non-uniform amplitude
  inside the crystal — a ramp or a dispersed front; a sharp step gives only the
  floor.

**Corrections to earlier claims.**
- Claims 22/23 and Entry 34: phrase the condition as Theta_1(final) = 0, not
  "S_f at a Bessel zero". Entry 34's tolerances (±5 % in S, 10 % loss -> 30:1)
  were for the MID case; for the fast ramp (T < t_tr) the tolerance is ~3x
  tighter (loss ~4 %).
- The step-onset control is not exactly zero during the transit: floor 0.0018
  at alpha / alpha_ac = 0.025.
- `nonideal.theta_law` uses c_{-n}: with a second harmonic c_{-n} != ±c_n, and the
  +n channel reads harmonic -n (as J_{-n} in Entry 14). A first run with c_{+n}
  gave law errors up to 0.3 at eps2 = 0.2; fixed and checked (error 0.014).

**Realism of the acoustic parameters (all E).**
- Attenuation of GHz shear waves in Si at room temperature, order 1–5 dB/cm
  (textbook order, not looked up) -> 0.3–1.7 % amplitude loss over 0.3 mm:
  inside the 30:1 region even for the fast ramp.
- Second harmonic eps2 ~ 0.01 (from beta ~ 5, strain ~ 7e-5, k L ~ 250) ->
  contrast ~2000 (FAST).
- Bulk-phonon dispersion at Lambda_ac ~ 5 um is negligible; the dispersion
  model stands in for any front broadening (transducer bandwidth, beam
  diffraction — not modelled).
- Ringing amplitude and decay are transducer-specific.

**Caveats.** All non-idealities are phenomenological; acoustic diffraction,
mode conversion and lateral inhomogeneity are absent. psi2 != 0 is checked by TT
only. Bragg was not repeated here (Entry 28 showed the same law; not re-checked
under non-idealities). Test:
`tests/test_formation.py::test_transient_only_contrast_set_by_accumulated_coupling`
(41 pass).

---

## 2026-10-01 — Entry 36: same amplitude, up vs down — the whole rocking curve (mostly a negative result)

**Question (user's review).** At two moments with the same transducer
amplitude, one on the way up a triangular pulse and one on the way down, the
depth profile inside the crystal differs. Does the *whole* I(eta) differ, and
how? Classification fixed in advance:
1. only relative satellite intensities differ -> history-dependent intensity;
2. widths / splitting / fringes differ visibly -> history-dependent line shape;
3. new peak positions or combination structure -> genuinely new.
Extra test: is each channel's line shape that of a *uniform* crystal with the
same accumulated coupling Theta_n (Laue slab of thickness 1, coupling
pi alpha Theta_n)? If yes, the depth history is compressed into one number per
channel.

**Setup.** `explore/nb/e36_updown.py`; figure `explore/make_updown_figure.py` ->
`figures/updown_rocking.png`. Laue, S_max = 3, alpha_ac = 40, eta step 0.25 over
±3.5 q, n_z from the extended rule. Per channel n = -3..3: window integral,
centre, FWHM, FWHM of the uniform-Theta line, an area-normalised shape
difference, local maxima; plus a search for peaks > 2 % of the maximum farther
than q/4 from any channel.

**A first run was confounded and is kept for the record** (`python e36_updown.py
partial`): with T = 1, A = 0.6 the "up" moment is t = 0.6, when the crystal is
only 60 % filled. Up-state widths of 11–16 (vs 5.5 down) and apparent "peaks
outside channels" at ±0.26–0.44 q (alpha = 2, T = 0.5) are the finite-depth sinc
fringes of a partly filled crystal (C-1.8), not history. **Rule: compare only
moments at which the crystal is fully filled (t_up = A T >= 1).**

**Result (both moments fully filled).** alpha = 2, T = 2, S A = 2.4
(t_up = 1.60, t_down = 2.40):

| n | window p up / down | centre up / down | FWHM up / down | uniform-Theta FWHM up / down |
|---|---|---|---|---|
| ±1 | 7.6–8.1 / 8.1–8.4 | 0.05–0.06 / 0.06 | 5.25 / 4.25–4.5 | 5.25 / 4.25–4.5 |
| ±2 | 8.0–8.1 / 4.9–5.2 | **0.99 / 0.08** | 5.75 / 4.25–4.5 | 5.0–5.25 / 4.5 |
| ±3 | 2.5–2.6 / 6.4 | **0.28 / 0.99** | 7.0–7.5 / 6.0 | 5.75 / 5.25 |
| 0 | 7.4 / 4.7 | 0.25 / 0.31 | 6.0 / 13.3 | 4.0 / 5.75 |

alpha = 2, S A = 1.8 (t_up = 1.20, t_down = 2.80): n = ±3 window 0.5–0.6 / 6.7,
centre 0.04 / 0.95; FWHM 8.5–9.25 / 6.0 (uniform 5.75 / 5.25–5.5).
alpha = 0.5, S A = 1.8: n = ±3 window 0.03–0.04 / 0.76; FWHM 8.5–9.5 / 6.25
(uniform 5.75 / 5.75); n = ±1 shape difference 0.05–0.07 (identical within noise).
Peaks outside channels: none, except one at 0.26 q in one case, a sinc fringe of
the broad n = 0 line (FWHM 11.5).

**Reading — classification.**
- **Level 1: yes, strongly.** Which satellite is brightest is different on the
  way up and down at the same transducer amplitude (S A = 2.4: n = ±2 on the
  way up, n = ±3 on the way down); integrated intensities differ up to ~11x
  (n = ±3, S A = 1.8).
- **Level 2: visible, but almost entirely compressible into Theta_n.** The n = 2
  line is single-peaked on the way up and split around a dark centre on the
  way down (figure b) — a striking difference, but the uniform-crystal line
  with the same Theta_2 reproduces both. What Theta_n does *not* capture is
  modest: up-state high orders are 20–60 % broader than their uniform-Theta
  line (n = ±3: 7.0–9.5 vs 5.75), and their side lobes are washed out. The
  reason is where the coupling sits: on the way up, the large amplitudes (which
  drive high orders) are only in the shallow, most recently filled layer, so the
  effective diffracting thickness is smaller and the line broader. Width
  therefore carries one extra piece of information — the depth extent of the
  coupling — beyond Theta_n.
- **Level 3: no.** No new peak positions, no combination structure.
- **Verdict:** "history-dependent satellite intensities (large, order-reordering),
  with mildly history-dependent line shapes." As the user anticipated, this is
  the resolved-channel structure at work: each channel is compressed to
  Theta_n plus a weak width correction. Not a reason to promote "memory" to the
  main story; it supports the travelling-front formation dynamics as the core.

**Possible follow-up (not done).** Inverse problem: the set {Theta_n} are Bessel
moments of the distribution of amplitudes inside the crystal,
Theta_n = int J_n(x) rho(x) dx; the widths add the depth extent. Whether a
pattern determines the recent drive history (and how uniquely — Bessel moments
do not fix the depth *order*) is a well-posed question for later.

**Caveats.** One pulse shape (triangle), S_max = 3, Laue only; FWHM on a 0.25 eta
grid (±0.25 resolution).

---

## 2026-10-01 — Entry 37: two acoustic domains — coherent cancellation at the centre, incoherent sum in the integral

**Hypothesis (user's proposal).** Fill the crystal with two acoustic domains of
amplitudes S1, S2 whose J_1 have opposite sign. Each alone diffracts strongly;
when their accumulated couplings balance, r J_1(S1) + (1 - r) J_1(S2) = 0, the
crystal is "dark" for n = 1. As the domain wall moves with the front, the
satellite should switch bright -> dark -> bright. Maximise the individual
strengths: S1 = 1.8412 (J_1 max, +0.5819), S2 = 5.3314 (J_1 min, -0.3461),
balance at wall depth r* = 0.3730.

**Setup.** `explore/nb/e37_domains.py`; figure `explore/make_domain_figure.py` ->
`figures/domain_cancellation.png`. Transducer: S2 for one transit, then S1, so for
t in [1, 2] the entrance-side layer zeta < t - 1 holds S1 and the rest S2.
"Domain alone" runs replace the other domain with an unmodulated crystal.
Laue and Bragg, alpha_ac = 40, alpha = 0.1 / 0.5 / 2, n_z from the extended rule.

**Result — centre.**

| geometry, alpha | S1 domain alone | S2 domain alone | both (wall at r*) | law max err n = 1, 2, 3 |
|---|---|---|---|---|
| Laue, 0.5 | 0.113 | 0.112 | 3.4e-6 | 0.008, 0.006, 0.003 |
| Laue, 2 | **0.959** | **0.960** | **8.1e-4** | 0.031, 0.027, 0.014 |
| Bragg, 0.5 | 0.109 | 0.106 | 1.8e-5 | 0.006, 0.005, 0.003 |
| Bragg, 2 | 0.771 | 0.767 | 8.1e-4 | 0.025, 0.021, 0.011 |

The minimum sits at t = 1.374 (predicted 1 + r* = 1.373). Dark width (I below 1 %
of the maximum) ~0.1 t_tr at alpha = 0.5; ~0.01–0.03 t_tr in Bragg at alpha = 2
(limited by the time grid). At alpha = 2 the Laue centre has further dark points
that are ordinary Pendellösung zeros (|Theta| alpha = 1), e.g. near zeta_b ~ 0.9 —
the cancellation is the one at r*.

**Result — the whole satellite.** At the cancellation (Laue, alpha = 2) the
centre is 8e-4 but the window integral is 17.6, *larger* than the S1 domain alone
(16.6); the line is a doublet with side maxima at eta - q = ±4.5, ±10.7, ...
Kinematic check (window integral minus the static tail):

| alpha | both | S1 alone | S2 alone | Parseval 2 pi (pi alpha)^2 int J_1^2 dzeta (both / S1 / S2) |
|---|---|---|---|---|
| 0.1 | 0.1226 | 0.0771 | 0.0463 | 0.1249 / 0.0783 / 0.0466 |
| 0.5 | 2.862 | 1.856 | 1.115 | 3.122 / 1.958 / 1.165 |

"Both" equals the sum of the two domains (0.1226 vs 0.1234). The same holds for
Entry 35's retuned-S* attenuation case: centre contrast 1.1e5 / 2.4e4 / 3.8e3
(10 / 20 / 40 % loss) but **integrated contrast only 27 / 7.9 / 2.2**.

**Order selection.** With (S1, S2) above, n = 1 cancels at r* = 0.373 while
Theta_2 = +0.077 and Theta_3 = +0.226 there (n = 3 bright: 0.121 at alpha = 0.5,
0.79 in Bragg alpha = 2). Other pairs: (1.0, 4.5) and (2.5, 4.2) cancel only n = 1;
(3.0, 6.5) has cancellation points for n = 1, 2, 3 at r* = 0.312, 0.387, 0.103 —
the wall position selects which order goes dark.

**Reading.**
- **Coherent cancellation is real and exact at the satellite centre**, in both
  geometries: two domains that each reflect 96 % (Laue) or 77 % (Bragg) at
  eta = q give < 1e-3 together. At exact resonance every layer rotates about the
  same axis, so the angles (accumulated couplings) add with sign.
- **It does not darken the satellite.** Away from eta = q the layer rotations do
  not commute and the cancellation fails. Kinematically the amplitude is the
  Fourier transform of the depth profile of the coupling,
  A(delta) ~ int J_1(S(zeta)) e^{i delta (1 - zeta)} dzeta, so the centre is its
  zeroth moment (Theta_1) while the eta-integral is, by Parseval,
  2 pi (pi alpha)^2 int J_1(S(zeta))^2 dzeta — a sum of squares that no
  arrangement can cancel. The two domains add *incoherently* in the integral.
  **The integrated satellite vanishes only where J_1 = 0 locally.**
- What the moving wall switches is the line shape: single peak -> centre-dark
  doublet -> single peak, with the integral roughly conserved. A detector with
  fine angular resolution at eta = q sees a sharp switch; an angle-integrating
  detector sees almost nothing.
- Order selection by wall position works in the same centre-only sense.

**Corrections to earlier entries.**
- Entry 35, section 6 ("the essence is the accumulated coupling, not the Bessel
  zero") holds **for the satellite centre only**. For the eta-integrated satellite
  the relevant quantity is int J_1^2 dzeta, so the transient-only satellite's
  integrated darkness at the end requires J_1(S_f) ~ 0 locally — the Bessel zero
  does matter there. Entry 35's integrated check (section 7) covered only the
  ideal and non-retuned cases; the retuned case was added here (27 / 7.9 / 2.2).
- The test docstring of `test_transient_only_contrast_set_by_accumulated_coupling`
  now says "centre". New test
  `test_two_domains_cancel_at_centre_but_add_in_the_integral` (42 pass).

**Caveats.** Ideal step between domains (a real transducer has a finite switching
time: the wall becomes a ramp, which leaves the centre cancellation intact but
changes the profile). No acoustic reflection at the wall (same material). S2 =
5.33 is a large drive (strain ~1e-4 for Si 440 at ~1 GHz, E).

---

## 2026-10-01 — Entry 38: chirped acoustic waves — a Landau–Zener plateau, no focusing

**Goal (user, overnight exploratory brief).** Look for transient diffraction
phenomena that are *not* centre-intensity changes, ordinary Pendellösung,
finite-depth sinc fringes, adiabatic Bessel-ladder scans, or Parseval-conserving
redistributions. Change the local acoustic wavevector instead of the amplitude.

**Framework.** `explore/nb/qmod.py`: `PhaseDrive(Drive)` with an arbitrary
acoustic phase theta(zeta, t) (phi = S A(s) sin theta, no core change) and
`kinematic()`, the first-Born line with the full e^{-i phi} (all harmonics) —
the control for "is it just a Fourier transform". Sanity: PhaseDrive with
uniform q equals Drive to 7e-15; kinematic equals TT to 0.04 % at alpha = 0.02.

**Two bounds that frame everything below (derived before running).**
In the kinematic limit the satellite amplitude is the Fourier transform of the
coupling profile c(zeta) = (e^{-i phi})_n. (i) The peak obeys
|int c dzeta|^2 <= (int |c| dzeta)^2: **a uniform phase maximises the peak** — no
acoustic profile can "focus" a satellite kinematically. (ii) The integral is
fixed by Parseval, 2 pi (pi alpha)^2 int |c|^2 dzeta, and does not depend on any
phase profile. So qualitatively new *spectral weight* can only be dynamical.

**A1 — spatial chirp** q(zeta) = q0 [1 + beta (zeta - 1/2)] (graded medium at fixed
frequency; the front is kept at zeta = t, an approximation), fully filled crystal
(t = 1.5), n = 1, S = 1.5 (J_1 = 0.558), alpha_ac = 40. `explore/nb/e38_chirp.py A1`.
Prediction: the local resonance sweeps at kappa = d eta_c / d zeta = q0 beta; the
satellite is a two-level Landau–Zener crossing with coupling a = pi alpha J_1, so a
flat plateau of height P = 1 - exp(-2 pi Gamma), **Gamma = a^2 / kappa**, over the
swept range (kinematic: 2 pi Gamma).

| alpha | beta | Gamma | LZ P | TT plateau | TT integral (uniform) | kinematic integral |
|---|---|---|---|---|---|---|
| 0.05 | 0 – 0.4 | — | — | — | **0.048 for every beta** | 0.048 |
| 0.5 | 0.2 | 0.015 | 0.092 | 0.096 | 4.58 (3.70) | 4.79 |
| 1 | 0.2 | 0.061 | 0.319 | 0.327 | 16.0 (6.62) | 19.2 |
| 1 | 0.4 | 0.031 | 0.175 | 0.181 | 17.4 (6.62) | 19.1 |
| 2 | 0.1 | 0.49 | 0.954 | 0.918 | 25.6 (10.1) | 76.7 |
| 2 | 0.2 | 0.24 | 0.785 | 0.754 | 40.0 (10.1) | 76.7 |
| 2 | 0.4 | 0.12 | 0.536 | 0.563 | 53.5 (10.1) | 76.6 |
| 4 | 0.1 | 1.96 | 1.000 | 0.987 | 33.1 (22.9) | 307 |
| 4 | 0.4 | 0.49 | 0.954 | 0.910 | 94.0 (22.9) | 306 |

Figure `figures/dephasing_extinction.png` (a: alpha = 2 lines; b: plateau vs Gamma).

- Kinematic limit: the chirp only spreads the line; the integral is conserved
  exactly (0.048 for all beta) — Parseval, as bound (ii) says.
- **Dynamical: the chirp raises the integrated satellite intensity** by up to
  5.3x (alpha = 2) and 4.1x (alpha = 4) over the uniform crystal, towards (never
  above) the kinematic value: the chirp removes extinction. The plateau height
  follows Landau–Zener with no fit for Gamma <~ 0.3. For Gamma >~ 0.5 the points
  scatter below the curve: the crossing is not complete inside the slab when
  q0 beta is not >> 2 pi alpha J_1 (e.g. alpha = 4, beta = 0.05: kappa = 12.6 vs
  2a = 14), and the line keeps Pendellösung structure.
- **Classification: DYNAMICAL** (kinematic model off by up to 9x in the
  integral). Dimensionless control: Gamma = (pi alpha J_n)^2 / (n q0 beta) =
  (sweep length of the resonance) / (satellite extinction length)^2 scale — the
  "detuning change over one extinction length vs Darwin width" the user asked
  for. Orders differ: Gamma_n ~ J_n^2 / n, so a chirp can be adiabatic for one
  order and diabatic for another.
- **Precedent (not verified in the literature tonight):** this is the
  Landau–Zener / Penning–Polder physics of X-rays in bent or strain-graded
  crystals (Takagi; Chukhovskii & Petrashen — names from memory). The acoustic
  version only replaces the bend by a chirp and the coupling by pi alpha J_n(S).
  Expect it to be known in kind.

**A2 — frequency chirp at the transducer** (theta = -Psi(t - zeta),
Psi' = q0 (1 + b s)): `e38_chirp.py A2`. The satellite centroid tracks the mean local
wavevector over the filled depth to ±1–3 in eta for b = 0.1, 0.3, alpha = 0.5 and 2
(e.g. alpha 0.5, b 0.3, t = 1.5: +75.7 measured vs +75.7 predicted); the line
broadens to the swept range once the crystal is full. **Kinematic, ordinary**: a
moving satellite that maps the chirp.

**E — focusing search** (theta = q0 (zeta - t) + c2 (zeta - 1/2)^2 + c3 (zeta - 1/2)^3,
13 x 9 grid, `e38_chirp.py E`): at alpha = 0.05, 0.5, 1 the best peak/integral and the
narrowest FWHM are at c2 = c3 = 0 — the uniform profile, exactly as bound (i) says.
At alpha = 2 a cubic phase raises peak/integral by 1.33x and narrows the envelope
from 16.0 to 11.9 — but only because the uniform alpha = 2 line at t = 1.5 is past its
first Pendellösung zero (split); its value (0.086) stays below what a uniform
crystal reaches at its own best thickness (0.153 at alpha = 1). **No focusing.**

---

## 2026-10-01 — Entry 39: acoustic domain walls — q-jump and phase jump

**Setup.** `explore/nb/e39_walls.py`. The transducer emits q1 (phase 0) for one
transit, then q2 (phase dphi); for t in [1, 2] the wall is at depth zeta_b = t - 1.
Controls: the **incoherent** sum of the two sides alone (other side unmodulated),
the **kinematic** model, and the wall **smoothed** over one and three acoustic
wavelengths. n = 1 unless stated, S = 1.5, alpha_ac = 40, static tail subtracted.

**B — q-wall at depth 0.5, q2 = q1 + 2 pi dc.**

| alpha | dc | npk both / incoh. / kin | integral both / incoh. | shape dev. TT vs kin |
|---|---|---|---|---|
| 0.05 | 0.1 – 10 | 1–4 / 1–2 / same as TT | 0.048 / 0.047–0.048 | <= 0.001 |
| 0.5 | 0.1 | 1 / 1 / 1 | 3.68 / 4.43 | 0.06 |
| 0.5 | 4 | 4 / 2 / 4 | 4.45 / 4.44 | 0.01 |
| 2 | 0.1 | 8 / 3 / 1 | **9.84 / 25.3** | 1.36 |
| 2 | 1 | 6 / 3 / 3 | 16.5 / 24.9 | 0.76 |
| 2 | 2 | 7 / 4 / 3 | 19.9 / 25.4 | 0.38 |
| 2 | 4 | 5 / 4 / 4 | 23.6 / 25.5 | 0.26 |
| 2 | 10 | 7 / 7 / 2 | 25.6 / 25.6 | 0.34 |

Smoothing over one wavelength changes nothing (e.g. alpha 2, dc 1: 16.48 vs 16.54).
- Kinematic: the line differs from the incoherent sum only by the coherent cross
  term between the two halves (centre doubled when dc -> 0, i.e. a thicker
  crystal); TT = kinematic. Ordinary Fourier optics.
- Dynamical (alpha = 2): when the two domains share a resonance (dc small) the
  coherent crystal is extinguished (9.8) far below two independent halves (25.3);
  detuning the domains by more than the satellite Darwin width (2 pi dc >~ 4 pi
  alpha_1 ~ 14) decouples them and the integral returns to the incoherent sum.
  **Same physics as A1** — dephasing lifts extinction — with a step instead of a
  gradient.
- Moving wall (`Bt`, dc = 1, alpha = 2): the integral rises from 9.8 (wall at an
  end) to 16.5 (wall in the middle) and back; the envelope width alternates
  15.9 / 5.4 / 15.5 / 5.4 / 15.9 as the effective sub-crystal thickness passes
  through the first Pendellösung maximum. Non-monotonic but explained by the
  sub-crystal thicknesses. No peak between the two satellite structures, no
  avoided crossing, no anomalously narrow feature.

**C — phase wall (same q), integrated satellite vs dphi, wall at depth 0.5.**

| alpha | n | dphi = 0 | pi/4 | pi/2 | 3pi/4 | pi | 3pi/2 | kinematic (all dphi) |
|---|---|---|---|---|---|---|---|---|
| 0.05 | 1 | 0.048 | 0.048 | 0.048 | 0.048 | 0.048 | 0.048 | 0.048 |
| 0.5 | 1 | 3.70 | 3.77 | 3.96 | 4.14 | 4.19 | 3.94 | 4.78–4.83 |
| 2 | 1 | 10.1 | 10.8 | 12.1 | 13.0 | 13.4 | 11.8 | 76.4–77.2 |
| 2 | 2 | 6.49 | 7.88 | **9.45** | 8.09 | **6.50** | 9.43 | 13.3–13.7 |
| 2 | 3 | 0.85 | 0.88 | 0.88 | 0.85 | 0.93 | 0.98 | 0.89–1.01 |

**Wall-depth scan, dphi = pi, n = 1** (`Cz`; figure panel c): at alpha = 2 the
integrated satellite rises to **2.58x** the no-wall value at zeta_b = 0.3 and 0.7,
with a dip to 1.33x at 0.5; at alpha = 0.5 at most 1.17x. The centre follows
|Theta_1| = alpha_1 |1 - 2 zeta_b| (Entry 30 law): maxima 0.97 at zeta_b ~ 0.28, 0.72,
zeros at 0.05, 0.5, 0.95.
Is it just two independent sub-crystals? Integral with the wall vs the sum of the
two sides alone: 23.2 / 20.0 (zeta_b 0.2), 26.1 / 23.3 (0.3), **13.4 / 25.7 (0.5)**,
26.1 / 23.2 (0.7). Not an incoherent sum: at zeta_b = 0.5 the two halves rotate in
opposite senses by equal angles and near resonance the second half undoes the
first (a Pendellösung "echo"), suppressing the satellite to half of two
independent halves; off the midpoint the wall adds 12–16 % over independence.
Smoothing: one wavelength leaves the effect (26.1 at zeta_b 0.3; 13.7 at 0.5);
three wavelengths weaken it (22.2; 11.8).

- **Odd / even:** the wall acts through n dphi. n = 2 is blind to a pi wall
  (6.50 vs 6.49) and maximally affected by pi/2 (1.46x); n = 1 is maximally
  affected by pi. n = 3 (alpha_3 = 0.12, kinematic) does not change.
- **Classification: DYNAMICAL** (kinematic integral independent of dphi). Needs
  alpha_n >~ 0.5.
- **Precedent (not verified tonight):** an acoustic phase wall is an "acoustic
  stacking fault / antiphase boundary" for the satellite with fault phase
  n dphi; stacking-fault contrast in dynamical X-ray diffraction (fringes vs
  fault depth, changed integrated intensity) is classic. Expect known in kind;
  the acoustic, order-selective, movable version may not be.
- No localised "defect-state" line was found: the wall produces doublets and
  fringe changes, not a new narrow feature at a new eta.

---

## 2026-10-01 — Entry 40: coarse automated anomaly search

**Setup.** `explore/nb/e40_search.py`, results `explore/nb/e40_search.json`,
figure `explore/make_search_figure.py` -> `figures/search_summary.png`.
250 cases = 25 profiles x alpha {0.05, 0.5, 1, 2, 4} x S {1, 2.5}, n = 1, filled
crystal (t = 1.5). Profiles: spatial chirp beta {0.02–0.4}; q-wall dc {0.3, 1, 2, 4}
at depth {0.25, 0.5, 0.75}; phase wall {pi/2, pi} at the same depths; amplitude
ramp T {0.3, 1}. Each case: full-TT and kinematic satellite lines (static
subtracted) reduced to features (peak count, envelope FWHM, peak, integral,
centroid, skew, centre) and scored against the uniform-q control at the same
(alpha, S): topology change, |ln integral ratio|, the dynamical part of it
(TT minus kinematic), narrowing, concentration (peak / integral).

**Result.**
- **Top 15 by dynamical change of spectral weight: all are chirps, q-walls or
  phase walls at alpha = 2–4**, integral 3.5–6.9x the uniform control while the
  kinematic integral changes by 0.000–0.036 (log). One mechanism (Entries 38–39).
- **Top by topology / narrowing / concentration:** (i) strong chirps (beta = 0.4)
  give 13–23 peaks and an 81-wide envelope — at every alpha including 0.05, i.e.
  Fourier spreading; (ii) "narrowing" cases (ramp T = 1 at alpha 4: FWHM 8.6 vs
  32.4; q-wall dc 4 at alpha 4: 6.1 vs 32.4; phase wall pi/2 at alpha 2: 4.5 vs 18)
  are all relative to a control that is itself split by Pendellösung at strong
  alpha; the profile reduces the effective accumulated coupling and returns a
  single line. Explained by Theta compression (Entry 36) / sub-crystal thickness.
- **Narrowest line in the whole search: 4.5** (phase wall pi/2, alpha 2) — above the
  dynamical-narrowing floor 0.71 x 5.566 = 3.96 of a uniform crystal (Entry 4).
  No super-narrow feature.
- **No profile reduced the integrated satellite by more than 16 %**
  (q-wall dc 1, alpha 4: 0.84x); phase profiles almost always increase it.
- No peak at a new eta position outside the swept local resonances; no
  avoided-crossing signature; no feature that strengthens with alpha other than
  the extinction relief.

---

## 2026-10-01 — Entry 41: validation, classification, and a conservative summary of Entries 38–40

**Validation** (`explore/nb/e41_validate.py`; integral / peak):

| candidate | base | 2 n_z | eta step / 2 | window x 1.5 | wall smoothed 1 wl / 3 wl | kinematic integral |
|---|---|---|---|---|---|---|
| uniform control alpha 2 / 4 | 10.11 / 22.91 | | | | | |
| chirp beta 0.4, alpha 2 | 54.30 / 0.726 | 54.30 | 54.30 | 56.79 | — | 77.3 |
| chirp beta 0.1, alpha 4 | 33.22 / 0.999 | 33.22 | 33.22 | 33.53 | — | 307 |
| pi wall zb 0.3, alpha 2 | 25.99 / 0.982 | 25.99 | 25.99 | 26.16 | 26.08 / 22.19 | 76.4 |
| pi wall zb 0.5, alpha 2 | 13.31 / 0.500 | 13.30 | 13.31 | 13.49 | 13.69 / 11.76 | 76.4 |
| q wall dc 2 zb 0.5, alpha 2 | 20.19 / 0.892 | 20.19 | 20.19 | 20.24 | 20.15 / 20.39 | 76.8 |

Converged to 4–5 digits in n_z and eta step; the window enlargement adds
< 1 % except for the widest chirp (+4.6 %, its plateau extends beyond the base
window). Smoothing the walls over one acoustic wavelength leaves every effect
intact; over three wavelengths the pi-wall gain drops from 2.6x to 2.2x (still
large). Not artefacts of sharp discontinuities. Kinematic integrals are the same
(~76.8 for alpha = 2) for every profile — Parseval — so every change of spectral
weight above is dynamical. Test:
`tests/test_formation.py::test_chirp_landau_zener_plateau_and_phase_wall_integral`
(43 pass).

**Kinematic vs dynamical classification.**

| effect | kinematic? | dynamical content | control parameter |
|---|---|---|---|
| chirp broadening, multi-peak envelope | yes (identical at alpha = 0.05) | none | q0 beta vs 5.566 / L |
| moving satellite under a frequency chirp | yes | none | b |
| q-wall interference between two halves | yes at alpha <~ 0.5 | none | dc |
| phase wall: centre-dark doublet | yes (Fourier) | — | n dphi |
| **chirp: integral x up to 5–7, Landau–Zener plateau** | no | **extinction relief by adiabatic passage** | Gamma = (pi alpha J_n)^2 / (n q0 beta) |
| **q-wall: integral recovers to the incoherent sum when 2 pi dc > Darwin width** | no | **decoupling of domains** | 2 pi dc / (4 pi alpha_n) |
| **pi phase wall: integral x 2.6 at zb 0.3, echo dip at zb 0.5** | no | **rotation reversal ("acoustic stacking fault")** | alpha_n, n dphi, zb |
| "focusing" / super-narrow line | no gain found | — | — |

**Physical interpretation (one sentence).** All new spectral weight found
tonight is one dynamical mechanism: a strong satellite (alpha_n >~ 0.5) is
extinction-limited, and anything that dephases the satellite resonance along the
depth — a gradient (chirp), a step in q (domain wall) or a step in phase (phase
wall) — breaks the coherent extinction and moves the integral towards its
kinematic (Parseval) value, by amounts set by Landau–Zener (gradients) or by
sub-crystal rotations (steps).

**Experimental plausibility (E unless stated).** A transducer frequency chirp
(A2) or a frequency step / phase step in an RF burst is standard; a spatial chirp
needs a graded crystal. For Gamma ~ 0.1–0.5 at alpha_1 ~ 1 (Si 440 at 17–30 keV,
L ~ 1–2 Lambda_P), q0 beta ~ 2 pi alpha_1 ... 20 alpha_1, i.e. a relative frequency
sweep of order 5–30 % across one transit — feasible for a broadband transducer
(E). The pi phase wall needs only an RF phase flip. All effects need resolved
channels (alpha / alpha_ac <~ 0.05, Entry 19) and alpha_n >~ 0.5.

**Strongest positive result.** The Landau–Zener plateau of a chirped satellite:
TT matches 1 - exp(-2 pi Gamma) without a fit over Gamma = 0.01–0.3 and four
values of alpha, with the integral raised up to 5.3x (alpha 2) — a quantitative,
dynamical, order-selective (Gamma_n ~ J_n^2 / n) effect.

**Strongest negative result.** No acoustic phase profile focuses a satellite:
kinematically impossible by Cauchy–Schwarz (proved), and dynamically no profile
gave a line narrower than the uniform crystal's own dynamical minimum (4.5 vs
floor 3.96) or a peak / integral above the uniform optimum. No transient peak at
a new eta, no avoided crossing, no defect-like line.

**What to test next.**
1. Literature: Landau–Zener / Penning–Polder in bent and graded crystals; stacking
   fault and antiphase-boundary contrast in dynamical X-ray theory — decide
   whether the acoustic, order-selective, time-dependent versions add anything
   publishable.
2. Transient versions: a chirped *front* entering the crystal (resonance-crossing
   point appearing inside the slab) — does the LZ plateau form in time t_tr / ...?
   (not run tonight beyond A2's centroid).
3. Bragg geometry for the chirp (bent-crystal Bragg reflectivity is known to rise
   with strain gradient; check the satellite version).

### Brutally conservative summary

**A. Genuinely interesting and worth pursuing:** none qualitatively new. The
closest candidate is the *order-selective, tunable* Landau–Zener extinction relief
of acoustic satellites (Gamma_n ~ (pi alpha J_n(S))^2 / (n q0 beta)) — quantitative and
dynamical, but its physics is very likely the known bent/graded-crystal
Landau–Zener problem. Worth one literature check before investing more.

**B. Real but probably ordinary:**
- chirp -> Landau–Zener plateau and integral increase (bent-crystal analogue);
- q-wall -> domain decoupling restores the incoherent sum (bilayer / mosaic analogue);
- pi phase wall -> integral x 2.6 and an echo dip at the midpoint, odd/even
  selectivity via n dphi (stacking-fault / antiphase-boundary analogue);
- moving satellite under a frequency chirp (maps the chirp; kinematic).

**C. Artefact / already explained:**
- multi-peak, wide lines under strong chirps (Fourier spreading, identical at
  alpha = 0.05);
- apparent narrowing and concentration at alpha >= 2 (control split by
  Pendellösung; profile reduces Theta);
- apparent focusing gain at alpha = 2 (same reason);
- centre-dark doublets from phase walls (Fourier; Entry 37's Parseval argument).

**D. Unresolved:**
- the pi-wall integral exceeding the sum of the two sides alone by 12–16 % off the
  midpoint (coherent cross-term; mechanism plausible, not derived);
- whether any of B survives in Bragg geometry and with absorption;
- the transient formation of the LZ plateau behind a chirped front.

**Nothing qualitatively new appeared tonight.** Every non-trivial effect reduces
to dynamical extinction being relieved by dephasing, with known analogues.

---

## 2026-10-02 — Entry 42: an exact Floquet solver for the uniform travelling wave

**Project (user brief, multi-hour run).** Map the breakdown of the resolved
Floquet-channel picture for the simplest field — uniform S, single q, perfect
sinusoid, no envelope tricks — by driving the satellite spacing into the
dynamical width. Crystal fully filled unless stated (steady state); Laue,
non-absorbing.

**Method.** `src/tratt/floquet.py`. Expand D0 = sum_k a_k e^{+ikq zeta},
Dh = sum_m b_m e^{-imq zeta}; the TT equations become a constant-coefficient
system x' = i H x with

    H[a_k,a_k] = -k q,   H[b_m,b_m] = eta + m q,
    H[a_k,b_m] = pi alpha J_{k+m}(S) e^{-i(k+m) q t},   H Hermitian.

The series is redundant at zeta = 0; a_k(0) = delta_k0, b_m(0) = 0 satisfies the
boundary conditions, so by uniqueness the truncated series converges to the exact
solution. Channel eta = n q = resonance of b_{-n} with a_0 (direct coupling
J_{-n}); couplings can be masked to remove individual paths. A cheaper solver for
scans is still full TT (RK4); the scan driver picks the cheaper of the two.

**Validation** (`tests/test_floquet.py`, 7 tests): Floquet vs full TT, max
abs difference 1.6e-6 – 1.3e-5 over (alpha, alpha_ac, S) = (1, 20, 1.5), (2, 5, 3),
(3, 2, 2), (5, 1, 4), (2, 0.5, 6) (the RK4 error level); K + 4 changes nothing
(<= 3e-7); |D0|^2 + |Dh|^2 = 1 to <= 3e-7; S = 0 reproduces the Laue formula to
1e-10; all couplings masked -> 0.

**Exact two-band structure** (`explore/nb/e46_bands.py`). The TT system is a 2x2
ODE with coefficients periodic in zeta (period 2 pi / q), so it has exactly two
Floquet exponents lambda_1, lambda_2 with lambda_1 + lambda_2 = eta (mod q),
obtained from the one-period monodromy (unitary: |det M| - 1 <= 6e-7). This is
the X-ray band structure of the acoustic superlattice; satellites are its gaps at
eta = n q. The truncated Floquet matrix carries copies of the same two branches
shifted by multiples of q.

---

## 2026-10-02 — Entry 43: the regime map — and most of the "breakdown" at one phase is micromotion

**Setup.** `explore/nb/e43_phase.py`, 1200 points: alpha in {0.05 ... 10} (10),
alpha_ac in {0.5 ... 40} (12), S in {0 ... 8} (10); full spectra over all orders
|n| <= ceil(S) + 4 (eta step 0.3); cached JSON per point in `cache_e43/` (acoustic
phase t = 0) and `cache_e43_ph4/` (average over 4 phases, `E43_NPHASE=4`).
Independent-channel prediction: every order a separate thickness-1 Laue slab of
coupling pi alpha |J_n| at n q, orders added in intensity. Metrics: E_pop
(order-window populations), E_spec (whole line), E_ctr (centres), participation
ratio PR (vs Bessel PR_ind), total diffracted intensity, peak count, flux
(max |I + I0 - 1| <= 1.7e-3 over all points; RK4 at the coarsest grids).
Figure `explore/make_phase_figure.py` -> `figures/regime_phase_diagram.png`.

**First (single-phase) map was misleading.** At t = 0, E_pop is 0.04–0.10 even
at R = 0.05 and 0.2–0.45 for alpha_ac < 2 at any alpha. Averaging over the acoustic
phase changes this completely:

| case | E_pop, t = 0 | E_pop, phase-averaged |
|---|---|---|
| alpha 0.25, alpha_ac 5, S 1.5 (R 0.05) | 0.043 | 0.001 |
| alpha 1, alpha_ac 5 (R 0.2) | 0.065 | 0.017 |
| alpha 2.5, alpha_ac 5 (R 0.5) | 0.217 | 0.145 |
| alpha 5, alpha_ac 5 (R 1) | 0.563 | 0.425 |

Reason: channel n carries the phase e^{i n q t}; cross terms between channels
(overlapping tails, or overlapping sincs when alpha_ac <~ 1) oscillate with the
acoustic phase and average out. **A phase-locked (stroboscopic) measurement sees
the coherent cross terms; a phase-averaging measurement restores additivity of
the orders.** The +-n asymmetry of single-phase spectra (up to 0.35 at R = 1) is
the same micromotion: it reverses sign over the acoustic period and averages to
~0 (checked at 8 phases), consistent with the third-order theory being +-n
symmetric (C_{-n} flips sign with J_{-n}).

**Phase-averaged result** (alpha_ac >= 2; median E_pop):

| R = alpha / alpha_ac = Lambda_ac / Lambda_P | <= 0.05 | 0.05–0.2 | 0.2–0.5 | 0.5–1 | 1–3 | > 3 |
|---|---|---|---|---|---|---|
| E_pop | 0.002 | 0.009 | 0.038 | 0.21 | 0.53 | 0.72 |
| total / independent | 1.00 | 1.00 | ~0.99 | 0.88 | 0.48 | 0.31 |

alpha_ac < 2 behaves the same after phase averaging (0.019 at R <= 0.05).
Regime boundaries (phase-averaged E_pop = 0.05 / 0.2 / 0.5): **R ~ 0.3 / 0.75 / 2**.
The boundaries move to larger R with S (single-phase medians: E_pop = 0.2 at
R = 0.57 for S <= 3, 0.95 for S >= 5), as expected from the natural variable
X = 4 R max_n |J_n(S)| (gap of the strongest order / half zone; Entry 46): binned
by X the per-S medians collapse (e.g. X 0.3–0.6: 0.05–0.09). Global fits:
log E_pop = 0.70 log R - 1.32 (scatter x1.7); adding alpha_ac and S: x1.5.
Phase-averaged, below the knee E_pop ~ 0.1 R.

**Proposed regimes** (phase-averaged):
- **A, resolved (R <~ 0.3):** independent orders, errors <~ 4 %.
- **B, perturbative overlap (0.3 <~ R <~ 0.75):** indirect couplings through
  forward-beam sidebands matter (Entry 44); dark orders reappear (Entry 45).
- **C, strong hybridisation (0.75 <~ R <~ 2):** gaps fold at q/2, lines leave the
  channel positions, dark windows (Entries 46–48); total intensity drops.
- **D, R >~ 2:** acoustic wavelength longer than the extinction length; the
  spectrum is no longer organised by orders (E_pop ~ 0.7).
Plus, for phase-locked measurements only, a kinematic coherence regime at
alpha_ac <~ 1–2 (orders overlap as sincs).

**Strong S** (Part 7; S up to 8 in the map, S = 12, 16 at alpha_ac = 5):
PR / PR_ind = 1.00 for R < 0.3 at every S, 0.8–1.3 at R > 1; at S = 12 / 16 the
rms order is 8.43 / 8.43 (R 0.2), 8.82 / 8.74 (R 1), 8.63 / 8.65 (R 2) and
11.2 / 11.2, 11.7 / 11.7, 11.7 / 11.6 vs the Bessel values — even where E_pop = 0.3.
Larger S is more robust (E_pop at R = 2: 0.30–0.35 vs the 0.53 median).
**Negative result: no localisation, diffusion or concentration in order space
beyond the Bessel distribution;** the breakdown is in line shapes and positions,
not in which orders are populated.

---

## 2026-10-02 — Entry 44: what fails first — a hierarchy of Floquet paths

**Setup.** `explore/nb/e44_mechanism.py`, `e44b_sidebands.py`
(`e44b_sidebands.json`); figure `figures/reduced_model_convergence.png`.
Reduced models in the Floquet basis vs the converged one (single phase t = 0):
IND (independent two-level channels, intensities added), STAR (a_0 + all b_m:
one shared forward beam, no forward sidebands), A1 / A2 (forward sidebands
|k| <= 1 / 2 added).

**Result (E_pop vs converged, S = 1.5).**

| alpha_ac | R | IND | STAR | A1 | A2 |
|---|---|---|---|---|---|
| 5 | 0.05 | 0.042 | 0.002 | 0.000 | 0.000 |
| 5 | 0.2 | 0.066 | 0.030 | 0.005 | 0.001 |
| 5 | 0.5 | 0.269 | 0.211 | 0.043 | 0.013 |
| 5 | 1 | 1.185 | 1.202 | 0.345 | 0.116 |
| 5 | 2 | 1.770 | 1.238 | 0.516 | 0.184 |

(alpha_ac = 2 and S = 4 rows follow the same pattern.) Total intensity:
IND and STAR over-predict by 1.2–2.5x at R >= 1; sidebands remove the excess.

Forward sidebands required for E_pop < 0.02 (all b_m kept):

| R | 0.02–0.1 | 0.2–0.35 | 0.5 | 0.75 | 1 | 1.5–3 |
|---|---|---|---|---|---|---|
| k_req, S = 1.5 | 0 | 1 | 2 | 3–4 | 4 | 4–5 (saturates) |
| k_req, S = 4 | 0 | 0–1 | 2–4 | 5 | 7 | > 8 |

**Reading — the order in which the picture fails.**
1. First (R <~ 0.1, single phase): the *coherence between off-resonant tails*
   of different orders, which the intensity sum drops. Kinematic, first order in
   R, and phase-dependent (it is what averages out in Entry 43). One shared
   forward beam (STAR) already fixes it.
2. Next (R ~ 0.2–0.5): **indirect paths through forward-beam Floquet sidebands
   a_{+-1}** — X-rays scattered into a neighbouring satellite and back into a
   sideband of the forward beam before reaching the channel. This is the same
   path that resurrects dark orders (Entry 45).
3. Then (R >~ 0.75): higher sidebands, roughly k_req ~ 4 R, saturating at the Bessel
   width of the drive (auto_K(S) ~ S + 3) — at that point every Floquet mode
   participates and no reduced model is useful.
The two-beam reduction per satellite (IND) is therefore not what fails first in
a phase-averaged measurement: it is exact to O(R) after averaging, and fails at
O(R^2) through the forward sidebands.
The three-mode model of Entry 14 (a_0, b_0, b_{-n} per channel) sits between IND
and STAR; it was not re-run here.

---

## 2026-10-02 — Entry 45: Bessel-zero resurrection — a third-order indirect coupling, and a renormalised zero

**Theory (derived before the runs).** H couples a's only to b's (bipartite), so
there is no second-order path from a_0 to b_{-n}. Löwdin elimination of the
off-resonant states (a_k at -kq, b_m at (m + n) q) gives at third order

    g_eff = pi alpha [ J_{-n}(S) - (pi alpha / q)^2 C_n(S) ],
    C_n(S) = sum_{k != 0, m != -n} J_{k-n} J_{k+m} J_m / (k (m + n)),

(pi alpha / q = R / 2); the remaining third-order terms are proportional to
J_{-n} and vanish at its zero. At S = j (J_n(j) = 0) the channel population is
sin^2|g_eff| ~ (pi alpha)^2 (R/2)^4 C_n^2, i.e. ~ alpha^6 / alpha_ac^4.
C_0(j_01) = +0.604, C_1(j_11) = -0.358, C_2(j_21) = +0.259.

**Result 1 — scaling, no fit** (`explore/nb/e45_zero.py scan`; figure
`explore/make_zero_figure.py` -> `figures/bessel_zero_resurrection.png`).
Channel population = max over eta near n q of |b_{-n}(1)|^2 (exact Floquet):

| n | alpha | alpha_ac = 40 | 10 | 5 | 3 | 2 | 1.5 | 1 |
|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 8.8e-8 (1.000) | 2.2e-5 (0.991) | 3.5e-4 (0.964) | 2.5e-3 (0.904) | 1.1e-2 (0.80) | 3.7e-2 (0.85) | 8.8e-2 (0.42) |
| 1 | 1 | 3.1e-8 (1.000) | 7.8e-6 (0.994) | 1.2e-4 (0.976) | 9.1e-4 (0.936) | 4.2e-3 (0.86) | 1.4e-2 (0.88) | 4.2e-2 (0.55) |
| 2 | 1 | 1.6e-8 (1.000) | 4.1e-6 (0.995) | 6.5e-5 (0.981) | 4.9e-4 (0.949) | 2.3e-3 (0.89) | 7.4e-3 (0.90) | 2.5e-2 (0.62) |

(ratio to the third-order prediction in brackets; alpha = 0.5 and 2 in the script
output follow the same pattern; alpha = 2, alpha_ac = 1.5 reaches 0.90 for n = 0.)
Exact to <1 % for R <= 0.1, within 10 % to R ~ 0.3, failing at R >~ 0.5 (higher
orders). **A nominally dark order becomes fully bright at R ~ 1.**

**Result 2 — the path** (`e45_zero.py paths`). Largest terms of C_n:
n = 1: (m, k) = (-2, -1), i.e. **a_0 -> b_{-2} (the n = 2 satellite) -> a_{-1} (first
forward sideband) -> b_{-1}**, weight J_2^2 J_3 = 0.068; n = 0: the symmetric pair
through the n = +-1 satellites; n = 2: through n = 3. Removal tests (alpha 1,
alpha_ac 3): removing that b_m or that a_k halves the population (the mirror path
remains); keeping only that path's three couplings gives 3.78e-5 vs its predicted
3.54e-5 (n = 1). The indirect coupling goes through the *outer* neighbour and the
forward sideband.

**Result 3 — renormalised zero** (`e45b_shifted_zero.py`). Because g_eff is real,
S alone can make it vanish, at S* = j + (R/2)^2 C_n(j) / J'_{-n}(j):

| n | R | S* full | S* predicted | population at S* (at j) |
|---|---|---|---|---|
| 0 | 0.20 | 2.39292 | 2.39319 | 1.4e-7 (3.5e-4) |
| 0 | 0.67 | 2.26814 | 2.27556 | 7.4e-5 (1.1e-1) |
| 0 | 1.00 | 2.06652 | 2.11397 | 2.0e-3 (3.0e-1) |
| 1 | 0.20 | 3.82263 | 3.82283 | 5.8e-8 (1.2e-4) |
| 1 | 0.67 | 3.72790 | 3.73305 | 6.1e-5 (4.7e-2) |

The darkest amplitude moves by ~R^2 and the darkness is restored 3–4 orders deeper
than at the bare zero.

**Result 4 — what an experiment sees** (`e45_zero.py phys` + inline check).
For integer alpha_ac the physical intensity at eta = n q equals the channel
population and the prediction (e.g. alpha 1: 4.93e-7 / 7.85e-6 / 1.24e-4 at
alpha_ac 20 / 10 / 5). For non-integer alpha_ac and in window integrals the
off-resonant tails (first order in R, phase dependent) dominate in the resolved
regime (window max 6e-4 vs 5e-7 at alpha_ac 20): the resurrection is a centre /
phase-locked observable there, and becomes a window-scale effect only at R >~ 0.5.
Gap picture (Entry 46): at J_1 = 0 the n = 1 gap reopens as 0.02, 0.34, 2.50, 5.87
for R = 0.2, 0.5, 1, 2.

---

## 2026-10-02 — Entry 46: Floquet bands, gap folding, mixing-point annihilation, peak trajectories

**Setup.** `explore/nb/e46_bands.py` (monodromy, two exact branches),
`e48_trajectories.py` (phase-averaged spectral maxima and the eta where a Bloch
mode is 50/50 forward / diffracted, alpha_ac = 2, S = 0.5, 1, 2.3, R = 0.05–2);
figure `explore/make_floquet_figures.py` -> `figures/floquet_structure.png`.

**Gaps at eta = n q vs the independent 2 pi alpha |J_n|** (alpha_ac = 5):

| S | R | n = 0 | n = 1 | n = 2 | n = 3 |
|---|---|---|---|---|---|
| 1.5 | 0.05 | 0.80 / 0.80 | 0.88 / 0.88 | 0.37 / 0.36 | 0.10 / 0.10 |
| 1.5 | 0.5 | 7.08 / 8.04 | 8.91 / 8.76 | 4.18 / 3.65 | 1.76 / 0.96 |
| 1.5 | 1 | 8.32 / 16.08 | 13.18 / 17.53 | 10.60 / 7.29 | 6.02 / 1.92 |
| 3.83 (J_1 = 0) | 0.2 / 0.5 / 1 | — | 0.02 / 0.34 / 2.50 | — | — |

**Reading.**
- The independent picture is "isolated gaps of width 2 pi alpha |J_n|". In the
  reduced zone a gap cannot exceed q/2; the measured gaps follow the independent
  line and **fold back at q/2** (n = 1 reaches 0.49 q at R ~ 0.85 and then
  decreases). The breakdown condition is therefore precisely
  **R |J_n(S)| ~ 1/2** for each order, i.e. X = 4 R max|J_n| ~ 2 for the spectrum —
  the natural variable of Entry 43. n = 0 departs earlier (max 0.28 q at R ~ 0.8)
  because more paths feed it.
- Weak orders gain gap through indirect paths (n = 3 at R = 1: 3.1x).
- **Mixing points:** at R -> 0 a Bloch mode becomes 50/50 only at eta = n q. With
  increasing R one mixing point moves out from eta = 0 (S = 1: 0.002, 0.048, 0.173,
  0.433 q at R = 0.05, 0.25, 0.45, 0.65) while another moves in from the zone edge;
  they meet and annihilate (S = 1 near R ~ 0.55; S = 0.5 around R ~ 1), after which
  no 50/50 point exists in that half-zone. (Points pinned at exactly q/2 are a
  branch-sorting artefact and are excluded.)
- **Peak trajectories:** phase-averaged maxima leave the channel positions:
  S = 1, the n = 1 peak moves 1.00 q -> 0.94 (R 0.45) -> 0.87 (0.65) -> 0.75 (0.85);
  the main peak at eta = 0 disappears for R ~ 0.85–1.25 (Entry 47). The spectrum
  reorganises around non-integer positions (e.g. S = 1, alpha 2, alpha_ac 2:
  maxima at +-0.617 q and +-1.373 q).
- Is it a true eigenmode hybridisation or overlapping finite peaks? The folding
  and the mixing-point annihilation are properties of the exact two-band
  structure, so the reorganisation is genuine band hybridisation, not a
  superposition artefact. It is, however, the generic behaviour of a two-band
  periodic (Hill-type) problem once gaps reach the zone size; nothing
  exceptional-point-like was found (the problem is Hermitian; branches never
  coalesce).

---

## 2026-10-02 — Entry 47: dark windows with strong coupling — multi-path interference at the X-ray–acoustic resonance

**Setup.** Candidates from the cache (darkest orders with |J_n| > 0.2);
`explore/nb/e47_dark.py` (diagnosis of alpha 2, alpha_ac 2, S 1, n = 0),
`e47b_darkmap.py` (maps of D_n = p_n / p_n^ind over R in [0.05, 1.5] x S in [0.2, 6],
alpha_ac = 2 and 5, n = 0, 1, 2, 4 phases), figure `explore/make_darkmap_figure.py`
-> `figures/darkmap.png`.

**Main-line dark state (alpha 2, alpha_ac 2, S 1; J_0 = 0.765, pi alpha J_0 = 4.8).**
Window |eta| < q/2 holds 0.0185 of the independent value (0.006–0.023 at every one of
8 acoustic phases). Sensitivity: sharp minimum in S (0.90: 0.164, 1.00: 0.019,
1.10: 0.022, 1.20: 0.131), alpha_ac (1.9: 0.153, 2.0: 0.019, 2.1: 0.130) and alpha
(2.2: 0.007, 2.5: 3.2). **Path removal:** no forward sidebands -> 1.18; no b_{+-1}
(the n = +-1 satellites) -> 1.05; no b_{+-2} -> 0.034 (still dark); a_0 + b_0 only ->
1.00. The darkness needs the direct path and the indirect paths through the
neighbouring satellites and forward sidebands, and disappears when either of them
is removed. Bands: the gap at eta = 0 is 6.23 ~ q/2 = 6.28 (maximal) and the two
Bloch modes there are 99.5 % / 0.5 % diffracted — unmixed, so the incident beam
projects on one mode and does not Pendellösung-oscillate. Where the weight went:
0.17 inside the window vs 4.79 in q/2 < |eta| < 3q/2 (independent total 19.1), with
maxima at +-0.617 q and +-1.373 q (Entry 46).

**Maps.** n = 0 dark valley at R ~ 1–1.35 for both alpha_ac (darkest D = 0.002 at
alpha_ac 5, R 1.05, S 0.8); its locus satisfies R J_0(S) ~ 0.8–1.2, i.e. the main-line
gap 2 pi alpha J_0 equals q — the main-line Pendellösung length equals the
acoustic wavelength. n = 1 dark valley at S ~ 2.3, R 1.15–1.35 (D 0.06; near the
J_0 zero), n = 2 at S ~ 3.4–3.6, R ~ 1.5 (D 0.05). Path removal for n = 1 (S 2.3,
R 1.25): removing the *outer* neighbour b_{-2} -> 0.46, the inner neighbour -> 0.14,
forward sidebands -> 0.99; n = 2 likewise (outer neighbour -> 0.65). The
third-order renormalised-zero branches (dotted in the figure) are separate from
these valleys. Bright counterparts: the Bessel-zero columns (D >> 1).

**Time domain** (Entry 49, case 5, alpha_ac 5, R 1.05, S 0.8): the main window
empties monotonically as the front fills the crystal: 21.4 (t = 0, unmodulated
reflection) -> 11.2 (0.5) -> 0.9 (0.9).

**Reading.** A genuine multi-path effect in the sense of the brief: individual
couplings are large, and the window darkens because the direct and indirect
paths interfere; removing one interfering path restores it. Physically it is the
regime where the acoustic wavevector matches the dispersion-surface gap of the
reflection — **the X-ray(neutron)–acoustic resonance k_s ~ Delta K_0** of Entin,
Iolin–Entin and Raitman et al. (Entry 51), here in a non-absorbing Laue crystal
with the gap renormalised by J_0(S). "Dark" means the reflection has moved out of
its angular window (lines displaced to non-integer positions), not that the
crystal stops diffracting.

---

## 2026-10-02 — Entry 48: time-resolved formation across the regimes

**Setup.** `explore/nb/e49_time.py`, travelling step front, alpha_ac = 5, window
integrals p_n(t) and centres, t in [0, 1) every 1/(16 alpha_ac); compared with the
independent formation law (window integral of a Laue slab of thickness t,
coupling pi alpha |J_n|); FFT of the traces. Figure `figures/time_regimes.png`.
Regimes: 1 resolved (alpha 0.5, S 1.5, R 0.1); 2 onset (2.5, 1.5, R 0.5); 3 strong
overlap (5, 1.5, R 1); 4 J_1 = 0 (5, 3.83, R 1); 5 dark main line (5.25, 0.8, R 1.05).

**A first version was wrong and is recorded:** it averaged over each acoustic
period (8 samples) to remove micromotion; the period is 0.2 t_tr, which
coarse-grained time to 0.2 and quantised every frequency and peak time.

**Result.** Relative deviation from the independent formation law (integrated
over the transit): regime 1: n = 1 0.075 (n = 0 2.3, n = 2 0.30 — weak/main windows
dominated by tails); regime 2: 0.19–0.78; regime 3: 0.33–0.58; regime 5: 0.31–0.77.
Spectral content (resolution ~1 cycle per t_tr — set by the transit itself):
regime 1: slow ~1 plus acoustic harmonics 5, 10, 15; regime 2: also 4, 6, 7, 11, 12
~ m alpha_ac +- alpha_n (alpha_n ~ 1.4); regime 3: 3 (~ alpha_n = 2.8), 7, 8, 14, 17 ~
m alpha_ac +- alpha_n.

**Reading.**
- In the overlap regime the formation clock alpha_n = R |J_n| alpha_ac approaches
  the acoustic clock alpha_ac; **slow formation and micromotion cease to be
  separable in time**, and traces carry Floquet combination frequencies
  m alpha_ac +- alpha_n. This is the time-domain image of the band structure
  (quasi-energy differences + m q), i.e. "ultrasound-induced Pendellösung
  beatings" (known in kind, Entry 51) — not a new effect.
- No non-Bessel formation order or population transfer between orders was
  found beyond what the band structure implies; the main-window emptying of
  regime 5 is the one qualitatively visible transient (Entry 47).
- Limitation: a single transit gives ~1 cycle per t_tr resolution; nearby
  frequencies cannot be separated in principle.

---

## 2026-10-02 — Entry 49: anomaly ranking against two reference models

**Setup.** `explore/nb/e50_anomaly.py` over the phase-averaged cache
(alpha_ac >= 2), results `e50_anomaly.json`, figure `figures/anomaly_summary.png`.
References: IND, and G3 = independent channels with the third-order couplings
|g_n| = pi alpha |J_{-n} - (R/2)^2 C_n| (level shifts and the J_{-n}-proportional
third-order terms omitted).

**Result.** Median E_pop, IND / G3: R <= 0.1: 0.0026 / 0.022; 0.1–0.3: 0.016 / 0.075;
0.3–0.6: 0.065 / 0.13; 0.6–1: 0.23 / 0.23; 1–2: 0.52 / 0.53; > 2: 0.66 / 0.81.
Top anomalies (largest E_g3, darkest strong orders, brightest near-zero orders) all
lie at R >= 1–1.5; the darkest strong orders are n = 0 windows (main-line
displacement, Entry 47: alpha 2, alpha_ac 2, S 1 is second) and the brightest
near-zero orders are J_2 ~ 0 at S = 5 with R 1.5–3.5, where windows are no longer
channels.

**Reading.**
- **Negative:** the naive third-order model is *worse* than IND away from Bessel
  zeros: the omitted level shifts and J_{-n}-proportional third-order terms are of
  the same order and compensate. G3 is valid at zeros only (where those terms
  vanish); a consistent global perturbative model was not built tonight.
- No anomaly at R < 1 survives beyond the mechanisms of Entries 43–47. Everything
  ranked high is regime C/D reorganisation.

(Entry 45 addendum: S* for n = 2: R 0.20: 5.12783 full vs 5.12799 predicted (population
3.0e-8 vs 6.5e-5 at j); R 0.33: 5.11342 vs 5.11442; R 0.67: 5.04666 vs 5.05082; R 1.00:
4.92908 vs 4.94483. n = 1, R 0.5: 3.77306 vs 3.77621; R 1.0: 3.58646 vs 3.60974.)

---

## 2026-10-02 — Entry 50: prior art for the overlap regime

Read in full: Polikarpov & Skadorov 1987 / 1988 (Entry 29); Raitman, Gavrilov &
Ekmanis, "Neutron diffraction on acoustic waves in a perfect and deformed single
crystals", InTech (2013), pp. 49–54 (theory section). Abstracts / snippets only:
Iolin & Entin, JETP 85, 1692 (1983); "Suppression of anomalous X-ray transmission
under the effect of periodic nonresonance deformations" (Crystallogr. Rep., ~2001);
acousto-optics (Klein–Cook parameter; deviations from the Raman–Nath Bessel law).

What is known:
- Raitman et al. state the regime structure explicitly: below the
  "neutron(X-ray)–acoustic resonance" (acoustic wavevector k_s below the
  dispersion-surface gap Delta K_0 = 2 pi / tau, tau = extinction length) the
  ultrasound is a perturbation with corrections ~ (H w)^2; near k_s ~ Delta K_0 the
  ultrasound mixes Bloch states, opens new gaps, causes **resonant ultrasonic
  suppression of the Borrmann transmission** and **ultrasound-induced Pendellösung
  beatings** (refs. to Entin, Iolin–Entin). Our R = Lambda_ac / Lambda_P ~ 1 is
  this resonance.
- Iolin–Entin and successors: superlattice period exceeding the extinction length
  (our R > 1) — analytical theory of Borrmann suppression in thick absorbing
  crystals.
- Polikarpov & Skadorov: multiwave Floquet treatment with J_n couplings, emphasis
  on non-resonant waves modifying extinction lengths; weak deformation xi << 1,
  three-wave truncation, steady state.
- Acousto-optics: finite Klein–Cook Q makes diffraction efficiencies deviate from
  the Raman–Nath Bessel law (zeros included) — a structural analogue of the
  resurrection, in a different (optical, thin/thick-grating) setting.

Classification of tonight's results:

| result | class |
|---|---|
| regime boundary at Lambda_ac ~ extinction length (R ~ 1) | **KNOWN** (X-ray/neutron–acoustic resonance) |
| suppression / emptying of the main reflection near resonance (Entry 47) | **DIRECT EXTENSION OF KNOWN** (non-absorbing Laue, J_0-renormalised condition 2 pi alpha J_0 ~ q, path-resolved) |
| combination frequencies in time traces (Entry 48) | **KNOWN** in kind (ultrasound-induced Pendellösung beatings) |
| gap folding at q/2, mixing-point annihilation (Entry 46) | **DIRECT EXTENSION** (standard two-band Floquet / Hill structure) |
| Bessel-zero resurrection with g_eff = pi alpha (R/2)^2 C_n, alpha^6/alpha_ac^4 law, dominant path through the outer neighbour and forward sideband (Entry 45) | **POSSIBLY NEW** for X-rays (Polikarpov's xi << 1 truncation cannot reach zeros; acousto-optic analogue exists) |
| renormalised Bessel zero S* - j = (R/2)^2 C_n / J'_{-n} (Entry 45) | **POSSIBLY NEW** (same caveat) |
| phase averaging restores order additivity; single-phase cross terms and +-n micromotion (Entry 43) | **UNCLEAR** (elementary, probably unremarked) |
| total diffracted intensity deficit collapsing on R (Entry 43) | **DIRECT EXTENSION** |
| no dynamical localisation in order space (Entry 43) | negative result |

Nothing here may be called new without reading Entin (JETP Lett. 26, 392 (1977);
JETP 77, 214 (1979)) and Iolin–Entin (1983) in full.

---

## 2026-10-02 — Entry 51: summary of the overlap-regime project

**What replaces the independent-satellite picture?** The exact two-band Floquet
structure of the X-ray in the acoustic superlattice. Satellites are its gaps; they
stay isolated and Bessel-weighted while R |J_n| << 1/2 (R = Lambda_ac / Lambda_P);
corrections come first from coherence between orders (phase-locked only), then
from indirect paths through forward-beam Floquet sidebands (O(R^2)); gaps fold at
q/2 when R |J_n| ~ 1/2, lines leave the channel positions, and near R ~ 1 the main
window empties at the X-ray–acoustic resonance. The population of orders keeps
the Bessel distribution throughout.

**A. Strongest genuinely interesting result.** Bessel-zero resurrection: a
nominally dark order is restored by a third-order indirect coupling
g_eff = pi alpha (R/2)^2 C_n(S), confirmed without a fit over four decades
(<1 % for R <= 0.1) for J_0, J_1, J_2 zeros, with the dominant path identified
(n = 1: through the n = 2 satellite and the first forward sideband) and a
renormalised zero S* - j = (R/2)^2 C_n / J'_{-n} (verified to <0.3 % at R <= 0.67).

**B. Strongest result likely already known.** The main-reflection dark window at
R J_0(S) ~ 1 (path-resolved multi-path interference, emptied as the front fills)
— the non-absorbing face of the known X-ray/neutron–acoustic resonance
(Entin; Iolin–Entin; Raitman et al.).

**C. Strongest negative result.** No dynamical localisation, diffusion or
concentration in Floquet-order space: the rms order equals the Bessel value to
two decimals even at S = 16, R = 2 where line shapes are badly broken. Also: the
naive third-order model is not a better global model (only at zeros), and the
single-phase "kinematic overlap" breakdown is micromotion that averages away.

**D. Cleanest scaling law.** Resurrected population P = sin^2[pi alpha (R/2)^2 C_n(S)]
~ alpha^6 / alpha_ac^4 (figure b of `bessel_zero_resurrection.png`). Second: total
diffracted intensity / independent prediction collapses on R for all S
(1.00 below R ~ 0.3, ~0.5 at R ~ 1–2). Regime boundaries (phase-averaged
E_pop = 0.05 / 0.2 / 0.5): R ~ 0.3 / 0.75 / 2, shifting to larger R with S as
X = 4 R max|J_n| suggests.

**E. Best candidate for the paper.** The validity map plus the resurrection law:
"the resolved Floquet picture holds for Lambda_ac < 0.3 Lambda_P (phase-averaged);
its first quantitative failure is the reappearance of Bessel-dark satellites with
a computable indirect coupling and a shifted dark amplitude". It fits the poster
as the boundary of the formation results (all earlier claims assumed
R <~ 0.05–0.1, now justified).

**F. Best candidate for a future project.** A consistent perturbative theory of
the overlap regime (including level shifts and the J-proportional third-order
terms) that reproduces the population errors for 0.1 < R < 0.5, and its comparison
with the Entin / Iolin–Entin resonance theory; plus the transient formation of a
resurrected channel (its clock is 1 / |g_eff|, i.e. very slow, so it lives on the
steady side of the transit).

**G. What should NOT be investigated further.** (i) Ranking dark/bright windows at
R > 1 — windows are no longer channels there; (ii) order-space localisation at
large S — absent; (iii) single-phase metrics as evidence of breakdown — they are
dominated by micromotion; (iv) caustic positions of outer peaks at large R by peak
picking (the simple S q + pi alpha law failed; the main-line Darwin width dominates).

**Deliverables.** Figures: `regime_phase_diagram.png` (1, 2), `bessel_zero_resurrection.png` (3),
`floquet_structure.png` (4), `darkmap.png` (5), `time_regimes.png` (6),
`anomaly_summary.png` (7), `reduced_model_convergence.png` (8). Code:
`src/tratt/floquet.py`; `explore/nb/e43_phase.py` ... `e50_anomaly.py`; caches
`cache_e43/`, `cache_e43_ph4/`. Tests: `tests/test_floquet.py` (7), 50 total pass.

---

## 2026-10-02 — Entry 52: effective 2x2 Hamiltonian of a satellite channel, derived independently and verified

**Starting point (exact).** Gauge t = 0. H[a_k,a_k] = -k q, H[b_m,b_m] = eta + m q,
H[a_k,b_m] = V J_{k+m}(S), V = pi alpha, H Hermitian. alpha = L / Lambda_P,
alpha_ac = L / Lambda_ac, q = 2 pi alpha_ac, R = alpha / alpha_ac = Lambda_ac / Lambda_P,
**V / q = R / 2** (the expansion parameter). Resonant subspace for channel n:
P = {a_0, b_-n} with E(a_0) = 0, E(b_-n) = delta = eta - n q; Q energies:
a_k = -k q (k != 0), b_m = (m + n) q + delta (m != -n). Independent-satellite limit:
H_IND = [[0, V J_{-n}], [V J_{-n}, delta]].

**Derivation** (quasi-degenerate perturbation theory, Schrieffer–Wolff / Winkler
form; denominators evaluated at E_P = 0, i.e. delta dropped inside them):
- Order 1: P-block of H.
- Order 2, diagonal: Delta_a = sum_{m != -n} |V J_m|^2 / (0 - (m + n) q) = -V (R/2) B_n;
  Delta_b = sum_{k != 0} |V J_{k-n}|^2 / (0 + k q) = +V (R/2) B_n, with
  **B_n = sum_{k != 0} J_{k-n}^2 / k** (substituting m' = m + n shows both sums are the
  same B_n). B_0 = 0 because J_{-k}^2 = J_k^2.
- Order 2, off-diagonal: **identically zero** — H is bipartite (a <-> b only), so no
  Q state couples to both a_0 (couples only to b's) and b_-n (couples only to a's).
- Order 3, off-diagonal, two kinds:
  (i) Q–Q paths a_0 -> b_m -> a_k -> b_-n: (1/2) sum V J_m . V J_{k+m} . V J_{k-n} x
  [1/((0 - E_bm)(0 - E_ak)) + same] = -V (R/2)^2 C_n,
  **C_n = sum_{k != 0, m != -n} J_{k-n} J_{k+m} J_m / (k (m + n))** (excluded states k = 0,
  m = -n);
  (ii) paths with one P-intermediate (wave-function renormalisation):
  -(1/2) V J_{-n} (R/2)^2 W_n, **W_n = sum_{m != -n} J_m^2 / (m + n)^2 + sum_{k != 0} J_{k-n}^2 / k^2**.
- Order 3, diagonal: identically zero (an odd path cannot return to the same
  sublattice).

    H_eff^(n) = [[ -V (R/2) B_n ,  g ], [ g , delta + V (R/2) B_n ]],
    g = pi alpha [ J_{-n}(S) (1 - (R/2)^2 W_n / 2) - (R/2)^2 C_n(S) ] + O(R^4),

so there are two distinct shifts: a **resonance-position shift in eta**,
delta* = Delta_a - Delta_b = **-pi alpha R B_n** (O(R^2) at fixed alpha_ac via alpha = R alpha_ac),
and a **coupling-zero shift in S**, g(S*) = 0 with
**S* - j = (R/2)^2 C_n(j) / J'_{-n}(j)** + O(R^4) (term (ii) does not enter at this order
because it is proportional to J_{-n}). At the bare zero g(j) = -pi alpha (R/2)^2 C_n(j),
so P_n ~ (pi alpha)^2 (R/2)^4 C_n^2 ~ alpha^6 / alpha_ac^4. The quoted Entry 45 result
is confirmed for the zero and the zero shift; for generic S it was missing term (ii).

**Numerical verification** (`explore/nb/e52_heff.py`; exact reference = half the
minimum splitting of the two Floquet branches nearest 0, found by bounded
minimisation over eta, and its eta position; alpha_ac = 10).

| n, S | R | rel. error g1 | rel. error 3rd (C only) | rel. error 3rd (full) | delta* exact / predicted |
|---|---|---|---|---|---|
| 0, 1.0 | 0.1 | 1.7e-3 | 7.6e-4 | 9.1e-6 | 0 / 0 |
| 1, 1.5 | 0.1 | 1.0e-4 | 5.6e-4 | **5.8e-8** | -0.11953 / -0.11962 |
| 1, 3.0 | 0.1 | 1.2e-3 | 3.3e-4 | 9.8e-7 | +0.01856 / +0.01847 |
| 2, 2.5 | 0.1 | 2.2e-5 | 3.8e-4 | 5.1e-9 | -0.10716 / -0.10712 |
| 2, 4.5 | 0.1 | 8.4e-4 | 1.7e-4 | 5.8e-8 | +0.01841 / +0.01839 |

Log-log slopes of the relative error over R = 0.01–0.1: first order and C-only
**2.00** in every case; full third order **3.98–4.09** (i.e. O(R^4), as derived).
delta* residuals scale with slopes 2.8–3.3 (next order), delta* itself ~R^2 at fixed
alpha_ac; n = 0 shows no shift, as B_0 = 0 requires.

At the Bessel zeros (n = 0, 1, 2; R = 0.02–0.2): |g| exact / predicted = 0.9999–1.0000
(R 0.02) down to 0.993–0.996 (R 0.2); **Delta S exact / predicted = 0.9994–1.0027**; at S*
the exact minimum splitting is 1e-11 – 7e-8, i.e. the gap closes completely — the
renormalised zero is a genuine level crossing, not a near-miss.

**Correction to Entry 49.** "The third-order model is worse than IND" came from two
errors: term (ii) was missing, and the G3 windows contained only each channel's own
line while IND sums all orders' lines (including their tails) before windowing.
Redone consistently (all orders summed, full third-order couplings, phase-averaged
cache, alpha_ac >= 2): median E_pop IND / full-third-order = 0.0026 / 0.0025 (R <= 0.1),
0.0162 / 0.0164 (0.1–0.3), 0.065 / 0.052 (0.3–0.6; improved in 78 % of cases). The
remaining population error at R <= 0.3 is not a coupling error (the coupling is now
exact to O(R^4)); it sits in the line shapes / tails, which the window model does
not describe.

Test: `tests/test_floquet.py::test_effective_2x2_hamiltonian_third_order` (51 pass).

---

## 2026-10-02 — Entry 53: is the acoustic phase average converged, and is the Bessel-zero residual visible in the exit intensity?

**Scope (user brief).** Verify, without new solvers: (1) whether the 4-phase average of
Entry 43 is converged; (2) whether the third-order residual coupling at a Bessel zero
(Entries 45, 52) dominates the *physical* exit intensity under some condition. External
numbers were treated as claims to check. No fifth-order theory, no new-phenomenon search.

**Definitions in the existing code.**
- Floquet: H[a_k,b_m] = pi alpha J_{k+m}(S) e^{-i(k+m) q t} (`src/tratt/floquet.py:63`);
  incidence x(0) = a_0 (`:82`); physical exit amplitude Dh(1) = sum_m b_m(1) e^{-i m q}
  (`:94`, zeta = 1); exit intensity |Dh(1)|^2.
- TT steady phase: phi(zeta) = S sin(q zeta - q t) (`src/tratt/drive.py`, `steady_phi`).
  Acoustic phase phi_ac = q t.
- Target population (Entry 45): `explore/nb/e45_zero.py:55–65`, |b_{-n}(1)|^2 at t = 0,
  maximised over eta in n q ± w. Channel eta = n q <-> b_{-n} (direct coupling J_{-n}).
- Entry 43 average: `explore/nb/e43_phase.py:62–64`, intensity averaged over
  t = k / (N alpha_ac), i.e. **phi_ac = 2 pi k / N, k = 0..N-1, N = 4**.
- Entry 45 window: |eta - n q| < 0.45 q (`e45_zero.py`, section `phys`).

**Derivation.** Let D = diag(e^{+i k phi_ac} on a_k, e^{-i m phi_ac} on b_m). Then
(D H(phi_ac) D^dagger)[a_k,b_m] = pi alpha J_{k+m} e^{-i(k+m)phi_ac} e^{+i k phi_ac} e^{+i m phi_ac}
= H(0)[a_k,b_m]; the diagonal is unchanged; D a_0 = a_0. Hence x(1; phi_ac) = D^dagger x(1; 0),
b_m(1; phi_ac) = e^{+i m phi_ac} b_m(1; 0), and

    A_h(phi_ac) = Dh(1) = sum_m c_m e^{i m phi_ac},   c_m = b_m(1; 0) e^{-i m q},  |c_m| = |b_m(1; 0)|.

So the uniform phase average is exactly <|A_h|^2> = sum_m |b_m(1)|^2 (Parseval over integer
m); **c_m equals the code's b_m up to the exit phase e^{-imq}, which drops out of |c_m|**.
An N-point uniform average keeps all cross terms with m - m' = 0 (mod N). Validity:
filled crystal, uniform travelling wave (steady state), non-absorbing, within the Floquet
truncation; *not* for the transient with a front inside the crystal. Checked numerically:
`FloquetLaue(t = phi/q)` reproduces e^{i m phi} b_m(0) to 1e-12, and the 16-phase TT
average equals sum |b_m|^2 to < 1e-3 (`tests/test_floquet.py::test_exact_acoustic_phase_average_is_sum_of_floquet_populations`).

**Reproduction** (`python explore/nb/e53_phaseavg.py check`; alpha = 1, S = j_11 = 3.8317,
eta = q; K = auto + 4 and K = auto + 14 identical):

| alpha_ac | target abs(b_-1)^2 | exact avg sum abs(b_m)^2 | N = 4 | N = 8 | N = 16, 32, 64 |
|---|---|---|---|---|---|
| 20 | 4.9226e-7 | 4.9226e-7 | 4.9234e-7 | 4.9225e-7 | 4.9226e-7 |
| 20.25 | 4.5939e-7 | 5.3856e-4 | 4.5268e-4 (abs -8.6e-5, -16 %) | 5.3919e-4 | 5.3856e-4 |
| 20.5 | 4.4632e-7 | 8.3845e-4 | 6.2890e-4 (abs -2.1e-4, -25 %) | 8.0193e-4 | 8.3845e-4 |

The external values are reproduced to all quoted digits. Direct TT at phi_ac = 0, pi/2, 3pi/2
equals the Floquet reconstruction to 4–5 digits (e.g. 8.5691e-4 / 8.5691e-4). Error vs N
(figure `figures/phaseavg_convergence.png`, `explore/make_phaseavg_figure.py`) falls
roughly exponentially (Bessel decay of c_m): ~1e-7 relative at N = 16, machine precision
for N >~ 22 (at this point, eta = q, S = j_11, alpha = 1 only). **The 4-phase average is not
converged at non-integer alpha_ac.**
*(Revised in Entry 54.)* At S = j_11, B_1(j_11) = -4e-17, so the **second-order resonance
shift** delta*_2 = -pi alpha R B_1 vanishes. This is not an exact coincidence of the resonance
with eta = q: the minimum two-branch gap lies at delta = -7.4e-6 (alpha_ac = 20), a resolved
higher-order shift (Entry 54). B_1 / J_1 varies with S (3.75, 0.68, -0.17, -0.07 at
S = 0.5, 1.5, 3, 5), so B_1 is not simply proportional to J_1 — why both vanish at j_11 is
not established.

**Impact on the Entry 43 map** (`e53b_map_check.py`, 16 representative points; then
`e53c_boundary.py` for the 277 points with 4-phase E_pop within ±30 % of 0.05 / 0.2 / 0.5,
and `e53d_affected.py` for all alpha_ac < 2 and all S >= 3: 769 points with 16 phases, cache
`cache_e43_ph16/`; merged E_pop in `e53_merged_epop.json`). Entry 43 does use a phase average
(4 phases); the reference is the independent-channel intensity sum, windows |eta - nq| < q/2.
- Representative points: cache reproduced exactly; exact average vs 4 phases differs by
  <= 0.009 (abs) at S = 1 but by 20–40 % (relative) at S = 5 (e.g. 0.031 -> 0.019,
  0.110 -> 0.085, 0.409 -> 0.328), always lowering E_pop.
- 769 points: **71 regime changes** (52 at alpha_ac < 2), 68 of them to a better regime;
  E_pop(16) - E_pop(4) median -0.005, range [-0.19, +0.06]. Largest: small alpha, alpha_ac 0.5–1,
  S >= 3: 4-phase 0.14–0.21 -> 16-phase 0.002–0.016.
- Medians (all alpha_ac >= 2 / alpha_ac < 2), 4-phase -> merged *(revised in Entry 54:
  scientific notation instead of values rounded to 0.000)*:
  R <= 0.05: 1.8e-3 / 1.9e-2 -> 1.7e-4 / 4.8e-4; 0.05–0.2: 8.8e-3 / 3.6e-2 -> 4.7e-3 / 4.2e-3;
  0.2–0.5: 3.8e-2 / 6.9e-2 -> 3.3e-2 / 3.9e-2; 0.5–1: 2.1e-1 / 1.9e-1 -> 1.9e-1 / 1.7e-1;
  1–3: 5.3e-1 / 4.7e-1 -> 5.0e-1 / 4.6e-1. Kinematic points (alpha <= 0.2): median 4.8e-3 -> 4.7e-4
  (max 2.1e-1 -> 7.4e-2).
- Boundaries (R where the binned median reaches 0.05 / 0.2 / 0.5; bins x1.3), 4-phase -> merged:
  S 0.5: 0.34/0.73/4.40 -> 0.44/0.73/4.40; S 1: 0.34/0.57/2.04 unchanged; S 2: 0.44/0.95/1.58 ->
  0.44/0.95/2.04; S 3 unchanged; **S 5: 0.20/0.95/2.04 -> 0.57/0.95/2.64; S 8: 0.20/1.22/2.64 ->
  0.57/1.22/2.64**. Points with S < 3 and alpha_ac >= 2 far from thresholds were not recomputed
  (representative checks there differ by <= 0.009).

**Bessel-zero residual in the exit intensity** (`e53_phaseavg.py eps`, alpha = 1, S = j_11,
alpha_ac = 20 + eps, -0.5 <= eps <= 0.5, refined to 1e-4 near 0; figure
`figures/bessel_residual_eps.png`). At eta = q (= corrected point here):
- target population ~4.5–5.5e-7, nearly independent of eps;
- exact phase-averaged exit intensity rises from 4.92e-7 (eps = 0) to 1.9e-4 (eps = ±0.1) and
  8–9e-4 (±0.5); the rest sum_{m != -1} |b_m|^2 grows ~eps^2 near 0;
- **target >= 50 % of the exit intensity only for |eps| <= 0.004, >= 90 % only for
  |eps| <= 0.0014** (contiguous intervals around the integer; on this eps grid — Entry 54's
  finer grid puts the edges at 0.00425–0.0045 and 0.00125–0.0015);
- window integral |eta - q| < 0.45 q: target fraction 4.4e-5 even at eps = 0 (3.1e-6 vs 7.0e-2):
  **not dominant in the integral over the examined window |eta - q| < 0.45 q** *(revised in
  Entry 54; narrower windows: Entry 54)*.
Integer-alpha_ac suppression is not exact: at eps = 0 the rest is 2.5e-13 (alpha 1, alpha_ac 20),
other/target = 3e-8 ... 1.3e-4 for alpha = 0.5–2, alpha_ac = 10, 20, scaling ~alpha^10 / alpha_ac^8
(i.e. ~alpha^2 R^8), unchanged with K = 30. Mechanism: the first-order off-resonant amplitudes
carry (e^{i(m+n)q} - 1), which vanishes for integer alpha_ac; higher orders remain. Shown for
S = j_11 only; not generalised.

**Verdicts on earlier claims.**

| claim | verdict |
|---|---|
| Entry 43: phase averaging restores order additivity | **maintained, sharpened**: exact average = sum_m abs(b_m)^2 |
| Entry 43: the 4-phase average is adequate | **revised**: not converged at non-integer alpha_ac or large S; use the exact formula sum_m abs(b_m)^2 (*revised in Entry 54*: 16 phases differ from it by up to 1.3e-2 in E_pop on the 437 checked points; no general N is claimed) |
| Entry 43: residual kinematic-overlap error at alpha_ac < 2 after averaging (0.019) | **revised**: 4-phase aliasing; merged average gives median 4.8e-4 (R <= 0.05; exact, Entry 54: 4.5e-4); the effect exists only for phase-locked measurements |
| Entry 43: regime boundaries R ~ 0.3 / 0.75 / 2 | **condition added**: S-dependent; converged values S = 0.5: 0.44 / 0.73 / 4.40; S = 1–3: 0.34–0.44 / 0.57–0.95 / 1.58–2.04; S = 5, 8: 0.57 / 0.95–1.22 / 2.64 (bin resolution x1.3; *corrected in Entry 54*: S = 0.5's 4.40 was omitted); not a universal law |
| Entry 45/52: third-order coupling and renormalised zero | **maintained** (target population reproduced exactly) |
| Entry 45: physical intensity at eta = nq equals the channel population at integer alpha_ac | **condition added**: true to other/target ~5e-7 at eps = 0 (alpha 1); dominance >= 50 % only within abs(eps) <= 0.004 and >= 90 % within 0.0014 (alpha 1, S = j_11); approximate, residual ~alpha^2 R^8 |
| Bessel-zero residual observable in angle-integrated intensity | **revised**: not dominant in the integral over abs(eta - q) < 0.45 q (fraction 4.4e-5); narrower windows: Entry 54 |
| B_1(j_11) = 0 | numerical observation (-4e-17), not proven; it removes the **second-order** shift only — the resonance itself is displaced by a resolved -7.4e-6 (Entry 54) |

**What the results support.** The phase average is exact as sum_m abs(b_m)^2 for the filled
uniform wave; the old 4-phase map was qualitatively right for alpha_ac >= 2 and small S but
overstated breakdown for small alpha_ac and large S. The third-order residual coupling is
the dominant exit intensity at the exact satellite centre only in a narrow window around
integer alpha_ac (relative tolerance of alpha_ac ~2e-4 for alpha = 1), and it is not dominant
in the integral over the examined window |eta - q| < 0.45 q *(revised in Entry 54; the
half-width dependence is in Entry 54)*.

**Remaining uncertainties / not done.** Angular resolution and thickness distribution of a real
experiment (no values given; left for a sensitivity analysis); other S and alpha for the eps
window; why B_1(j_11) = 0; the not-recomputed S < 3, alpha_ac >= 2 points far from thresholds.

**Reproduction.** `python explore/nb/e53_phaseavg.py check eps`;
`python explore/nb/e53b_map_check.py`; `python explore/nb/e53c_boundary.py`;
`python explore/nb/e53d_affected.py`; `python explore/make_phaseavg_figure.py`.

---

## 2026-10-02 — Entry 54: exact phase-averaged map, resonance position at the Bessel zero, target share vs integration half-width

**Scope (user brief).** Follow-up to Entry 53. It reuses the existing solvers and caches.
1. Restrict Entry 53's conclusions to the computed range (edited in place, marked "revised in
   Entry 54").
2. Check the phase-averaged map against the exact average sum_m |b_m|^2.
3. Locate the actual resonance at S = j_11.
4. Compute the target share of the exit intensity vs the angular integration half-width w and
   the period mismatch eps.

This is a signal-share analysis for one condition (alpha = 1, S = j_11, n = 1). It does not
establish experimental detectability or novelty. No experimental angular resolution or
thickness was assumed.

**1. Corrections to Entry 53 (in place).**
- "never dominant in an angle-integrated signal" -> "not dominant in the integral over the
  examined window |eta - q| < 0.45 q".
- B_1(j_11) ~ 0 is now stated as: the **second-order** resonance shift vanishes. This is
  separate from the resonance position itself (section 3).
- Regime-boundary verdict: added S = 0.5 (0.44 / 0.73 / 4.40). The value was in the bullet
  list but missing from the summary row.
- Medians are now in scientific notation (e.g. R <= 0.05: 1.7e-4 / 4.8e-4, not 0.000).
- Entry 53's point-dominance edges (|eps| <= 0.004 / 0.0014) were limited by the eps grid. The
  finer grid of section 4 puts them at 0.00425–0.0045 (50 %) and 0.00125–0.0015 (90 %).

**2. Exact phase-averaged map.**
- Code: `explore/nb/e54_exactmap.py`, cache `cache_e43_exact/`, output `e54_exactmap.json`;
  summary `explore/nb/e54c_summary.py`.
- Method: same eta grid, windows and independent-channel reference as `e43_phase.metrics`.
  The intensity is sum_m |b_m(eta)|^2 from `FloquetLaue.propagate`.
- Points: all S = 5 and 8 (240) and all points whose merged E_pop is within ±30 % of
  0.05 / 0.2 / 0.5. Total 437. All 437 also have a 16-phase value.

Results:
- Truncation: K and K + 4 differ by at most 1.0e-8 in E_pop.
- Exact minus 16-phase: median -3.0e-5; 99th percentile of the absolute value 2.9e-3;
  maximum absolute value 1.3e-2.
- The largest differences are at small alpha_ac and large S, (alpha, alpha_ac, S):
  (0.05, 0.5, 8): 1.5e-2 -> 1.8e-3; (0.1, 0.5, 8): 1.6e-2 -> 7.0e-3;
  (0.2, 0.5, 8): 3.2e-2 -> 2.8e-2.
  At these points the 16-phase average still has order-1 relative error in a small E_pop
  (aliasing of the m - m' = 16 terms). Elsewhere differences are <~ 3e-3.
- Regime changes, 16-phase -> exact: 2 of 437. Both straddle 0.05 by < 2e-4:
  (1, 2, 4): 5.01e-2 -> 4.996e-2; (5, 12, 5): 5.01e-2 -> 4.994e-2.
- Boundaries (bins x1.3) and medians do not change at the quoted precision.

| S | R at E_pop = 0.05 / 0.2 / 0.5: 4-phase | 16-phase merged | exact where available |
|---|---|---|---|
| 0.5 | 0.34 / 0.73 / 4.40 | 0.44 / 0.73 / 4.40 | 0.44 / 0.73 / 4.40 |
| 1 | 0.34 / 0.57 / 2.04 | 0.34 / 0.57 / 2.04 | 0.34 / 0.57 / 2.04 |
| 2 | 0.44 / 0.95 / 1.58 | 0.44 / 0.95 / 2.04 | 0.44 / 0.95 / 2.04 |
| 3 | 0.34 / 0.95 / 1.58 | 0.34 / 0.95 / 1.58 | 0.34 / 0.95 / 1.58 |
| 5 | 0.20 / 0.95 / 2.04 | 0.57 / 0.95 / 2.64 | 0.57 / 0.95 / 2.64 |
| 8 | 0.20 / 1.22 / 2.64 | 0.57 / 1.22 / 2.64 | 0.57 / 1.22 / 2.64 |

- Medians (alpha_ac >= 2 / alpha_ac < 2), exact where available: R <= 0.05: 1.4e-4 / 4.5e-4;
  0.05–0.2: 4.7e-3 / 4.1e-3; 0.2–0.5: 3.3e-2 / 3.9e-2; 0.5–1: 1.9e-1 / 1.6e-1;
  1–3: 5.0e-1 / 4.6e-1. Kinematic points (alpha <= 0.2): median 4.4e-4, maximum 7.4e-2.
- Not recomputed exactly: points far from thresholds with S < 5. They keep their 16-phase or
  4-phase values. Whether 16 phases suffice depends on the point. No general N is claimed.

**3. Resonance position at S = j_11, n = 1, alpha = 1.**
- Code: `e54b_width.py res`.
- Method: same as `e52_heff.exact`. Minimise the half-splitting of the two eigenvalues of
  `fl.H` nearest zero over eta = q + delta, with a bounded search in ±0.05.

| alpha_ac | delta* (min gap) | delta*/q | half-gap g(delta*) | g(0) - g(delta*) |
|---|---|---|---|---|
| 20 | -7.38e-6 | -5.9e-8 | 7.01928175e-4 | 9.70e-9 |
| 20.25 | -7.11e-6 | -5.6e-8 | 6.84708663e-4 | 9.23e-9 |
| 10 | -5.917e-5 | -9.4e-7 | 2.80517483e-3 | 1.55e-7 |

- Search and truncation error: K = auto + 4 / + 8 / + 12 and xatol 1e-9 / 1e-12 change
  delta* by <= 1.2e-8 (alpha_ac = 20) and g by <= 2e-13.
- Consistency: for two levels, g(0) - g(delta*) = delta*^2 / (8 g) = 9.7e-9 at alpha_ac = 20.
  This matches the direct difference. The excess of 9.7e-9 is far above the 2e-13 noise.
- **The difference is resolved. The resonance is not at eta = q.**
- The second-order shift -pi alpha R B_1 vanishes (B_1 = -4e-17), but a higher-order shift
  remains. At fixed alpha, the shift grows x8 from alpha_ac = 20 to 10. This is consistent
  with ~R^3 (fourth order). It is not derived here; no fourth- or fifth-order theory was
  attempted.
- Size: |delta*| = 7.4e-6, i.e. 1e-2 of the half-gap g and 6e-8 of q.

**4. Target share vs integration half-width.**
- Code: `e54b_width.py base eps map`.
- Figure: `figures/target_fraction_map.png` (`explore/make_phaseavg_figure.py`).

Definitions:
- Centre c: nominal (eta = q) or numerical resonance (q + delta*).
- Target I_t(w): integral of |b_{-1}|^2 over |eta - c| < w.
- Total I_a(w): integral of sum_m |b_m|^2 over the same window. This is the exact
  phase-averaged physical exit intensity (Entry 53).
- Share F(w) = I_t / I_a.
- Grid: offsets u = 0 plus geomspace(1e-8 q, 0.45 q, npts) on both sides of c; trapezoid rule.
- w_50, w_90: largest w for which F >= 0.5 (0.9) holds contiguously from the centre.

alpha_ac = 20 (eps = 0), q = 125.66:

| w/q | 1e-6 | 1e-5 | 1e-4 | 1e-3 | 3e-3 | 1e-2 | 3e-2 | 1e-1 | 0.45 |
|---|---|---|---|---|---|---|---|---|---|
| F | 1.0000 | 0.9997 | 0.974 | 0.270 | 3.9e-2 | 3.8e-3 | 6e-4 | 2e-4 | 4.4e-5 |

- **F >= 90 % for w/q <= 2.0e-4** (w = 0.025). There I_t = 2.5e-8 and I_a = 2.8e-8.
- **F >= 50 % for w/q <= 6.0e-4** (w = 0.076). There I_t = 7.4e-8 and I_a = 1.5e-7.
- For comparison, I_a(0.45 q) = 7.0e-2. The window where the target dominates holds about
  2e-6 of the exit intensity integrated over the Entry 45 window.
- Convergence: npts 1500 -> 3000 changes w_50 and w_90 by <= 0.3 %. The log-grid step
  (~1.2 %) limits resolution. K + 4 gives identical values.
- Nominal vs numerical centre: identical to 4 digits, because |delta*| = 6e-8 q is far below
  every w examined.

Period mismatch: centre at the numerical resonance for each alpha_ac; eps step 2.5e-4 in
|eps| <= 0.008, plus a log grid out to 0.05.

| abs(eps) | F(w -> 0) | w_90 / q | I_a(w_90) | w_50 / q | I_a(w_50) |
|---|---|---|---|---|---|
| 0 | 1.000 | 2.0e-4 | 2.8e-8 | 6.0e-4 | 1.5e-7 |
| 0.0005 | 0.987 | 1.9e-4 | 2.6e-8 | 6.0e-4 | 1.5e-7 |
| 0.001 | 0.952 | 1.5e-4 | 2.0e-8 | 5.9e-4 | 1.4e-7 |
| 0.00125 | 0.927 | 1.1e-4 | 1.5e-8 | 5.8e-4 | 1.4e-7 |
| 0.0015 | 0.898 | none | – | 5.7e-4 | 1.4e-7 |
| 0.002 | 0.833 | none | – | 5.4e-4 | 1.3e-7 |
| 0.003 | 0.689 | none | – | 4.5e-4 | 1.1e-7 |
| 0.004 | 0.554 | none | – | 2.7e-4 | 6.6e-8 |
| 0.00425 | 0.524 | none | – | 1.8e-4 | 4.5e-8 |
| 0.0045 | 0.496 | none | – | none | – |

- For |eps| >= 0.0015 no w gives F >= 90 %. For |eps| >= 0.0045 no w gives F >= 50 %. In
  both cases even the w -> 0 share is below the threshold.
- +eps and -eps agree to the printed digits (largest difference: F(w -> 0) 0.526 vs 0.524 at
  0.00425).
- F(0.45 q) = 4.4e-5 for all |eps| <= 0.05.
- delta* stays within -7.33e-6 ... -7.44e-6 for |eps| <= 0.05.

**Interpretation (this condition only).**
- The third-order residual dominates the integrated exit intensity only in a window
  w/q <~ 2e-4 (90 %) or 6e-4 (50 %) around the resonance.
- This holds only for |eps| <~ 1.3e-3 (90 %) or 4.3e-3 (50 %).
- The integrated intensity inside that window is ~1e-7 in `e43_phase` units, about 2e-6 of
  the Entry 45 window integral.
- I_t and I_a are integrals over the dimensionless eta. They are not photon efficiencies, and the
  2e-6 ratio is not a photon-loss rate. Converting them needs the incident angular distribution
  and the detector acceptance, which were not given.
- Measurability depends on angular resolution, flux and thickness spread. None were given;
  none were assumed.

**Verdicts.**

| claim | verdict |
|---|---|
| Entry 53: 16-phase map adequate | **condition added**: within 1.3e-2 of exact on the 437 checked points; 2 regime changes at the 0.05 edge; boundaries and medians unchanged; not a general N criterion |
| Entry 53: corrected resonance = nominal at S = j_11 | **revised**: only the second-order shift vanishes; resonance at delta* = -7.4e-6 (resolved, higher order) |
| Entry 53: residual never dominant in an angle-integrated signal | **revised**: not dominant for w = 0.45 q; >= 90 % for w/q <= 2.0e-4 and >= 50 % for w/q <= 6.0e-4 at eps = 0 |
| Entry 53: point-dominance window abs(eps) <= 0.004 / 0.0014 | **refined**: edges at 0.00425–0.0045 / 0.00125–0.0015 |

**Remaining / not done.** Other S, alpha and n; the origin of the R^3-like shift (not
derived); experimental resolution, flux and thickness distribution; detectability; novelty.

**Reproduction.** `python explore/nb/e54_exactmap.py` (parallel, ~25 min);
`python explore/nb/e54c_summary.py`; `python explore/nb/e54b_width.py base eps map res`;
`python explore/make_phaseavg_figure.py`.

---

## 2026-10-02 — Entry 55: literature check of the Entry 39 acoustic walls against stacking-fault and deformed-crystal dynamical theory

**Scope (user brief).** Before any new scan, decide whether the Entry 39 wall effects are
(a) a known stacking-fault or bicrystal effect, reproduced with an acoustic wave, or
(b) a new measurement or control capability, because the boundary moves and is tunable.
No calculations were run.

**Evidence level.** *(Superseded in part by Entry 56: Authier & Simon 1968 and Manukyan 2014 were later read in full and compared formula by formula.)* Search results, abstracts and bibliographic records only. The full texts
of Authier & Simon (1968), Kato, Usami & Katagawa (1967), Authier & Patel (1975) and Manukyan
(2014) were **not read**: they are paywalled or the fetch was refused (ITC Vol. B ch. 5.1:
HTTP 403). Every "known" verdict below rests on the stated scope of these works. None rests on
a formula-by-formula comparison.

**Sources and the scope they state.**
- A. Authier & D. Simon, Acta Cryst. A 24, 517–526 (1968).
  - Applies Takagi's dynamical theory to a single stacking fault in a perfect crystal.
  - Covers both plane-wave and spherical-wave incidence; fringe contrast via the stationary
    phase method.
  - The plane-wave theory for faults was later reformulated for the general Laue case.
- N. Kato, K. Usami, T. Katagawa, Adv. X-Ray Anal. 10 (1967), "The X-ray diffraction image of
  a stacking fault": spherical-wave Pendellösung theory extended to a crystal with one stacking
  fault in an arbitrary position.
- A. Authier & J. R. Patel, phys. stat. sol. (a) 27 (1975): X-ray topographic determination of
  whether stacking faults are intrinsic or extrinsic. The method relies on the sign of the
  fault phase.
- H. M. Manukyan, Izv. NAN Armenii Fizika 49, 360–365 (2014). This is the closest analogue
  found.
  - Dynamical X-ray diffraction by a **superlattice with a stacking fault**.
  - Gives the reflectance "within satellite limits" as a function of the wave's **phase jump**
    at the fault and the **fault depth**.
  - Whether the phase jump is resolved per satellite order is not stated in the abstract.
- Background: Bonse, Lang, Kato, Patel and Taupin developed dynamical diffraction in deformed
  crystals, and Takagi–Taupin theory is standard for faults and strain fields.
- Background: Polikarpov (pss(b) 143, 1987) gives a dynamical theory with ultrasound, steady
  state. Zolotoyabko controlled rocking curves with standing acoustic waves (quartz,
  15–105 MHz). Neither reference, as found, treats a phase or q discontinuity of the drive.

**Mapping.** Satellite n sees the coupling pi alpha J_n(S) e^{-i n phi_ac}. An acoustic phase
jump dphi at depth zeta_b is therefore, for satellite n, a planar fault with fault phase
alpha_f = n dphi at depth zeta_b. This is the same role as alpha_f = 2 pi h.R for a stacking
fault with displacement R. A q-jump changes the satellite's resonance position on one side,
which is the analogue of a lattice-parameter step (bicrystal or heterostructure).

| Entry 39 result | standard counterpart | verdict |
|---|---|---|
| C: integrated satellite depends on dphi, but only dynamically (kinematic integral independent of dphi) | stacking-fault contrast vanishes kinematically for a single fault in an integrated measurement; dynamical fault contrast is classic (Authier & Simon; Kato et al.) | **known in kind** |
| C: n = 2 blind to a pi wall, n = 1 most affected by pi; "acts through n dphi" | fault invisibility when h.R is an integer (alpha_f = 2 pi m); alpha_f depends on the reflection | **known in kind** (n plays the role of h) |
| C: fault-depth dependence; pi wall at mid-depth suppresses the satellite ("echo", Theta = alpha_1 abs(1 - 2 zeta_b)); x2.58 at zeta_b = 0.3 / 0.7 | plane-wave fault theory at arbitrary depth and phase (Authier & Simon; general Laue reformulation). The mid-depth pi case is a special case: coupling sign flips, and the second half undoes the first | **known in kind**; the specific numbers follow from the standard formula (not checked against it) |
| B: q-wall; shared resonance gives extinction, detuning > Darwin width decouples, integral returns to the incoherent sum | two crystals or layers with different lattice parameter in series; decoupling once the detuning exceeds the reflection width; extinction relief by dephasing (Entry 38 analogue: bent or graded crystals) | **known in kind** |
| B: moving wall, integral non-monotonic through sub-crystal Pendellösung | fault-depth fringes (the fault position sets the sub-crystal thicknesses) | **known in kind** |
| No new localised line at the wall | consistent with fault theory (fringes and doublets, no bound state for a single phase fault in two-beam theory) | consistent |
| Superlattice satellites with a fault phase jump | Manukyan (2014), static superlattice | **direct precedent** in a static system |

**What is not covered by the sources found (candidate differences; not established).**
1. *Tunable fault phase.* An acoustic fault phase alpha_f = n dphi is set electronically and
   can be varied continuously in situ. A crystallographic fault has a fixed R.
2. *Moving fault.* An acoustic fault moves at the sound speed, so the fault depth is a time
   variable. One sample then scans zeta_b stroboscopically. In topography, depth dependence
   comes from wedge-shaped samples or the fault inclination, i.e. a spatial scan.
3. *One dphi, all orders.* alpha_f = n dphi for every satellite at once, the analogue of
   different h for one R. Its use as a control knob (e.g. the n = 2 blind / n = 1 maximal
   pair) is not described in the sources found. In kind it is the h-dependence of fault
   contrast.
4. *A possible measurement use (idea only).* With alpha_f controlled, the fault contrast
   depends on the relative sign and phase of the couplings on the two sides. This could in
   principle read the sign of J_n(S) or of the third-order coupling (Entries 45, 52), which
   rocking-curve intensities do not show. Not checked; it may already follow from standard
   fault theory with complex structure factors.

**Verdict.** The physics of Entry 39 is the established theory of dynamical diffraction by a
planar fault or a lattice-parameter step, with the acoustic phase or q as the fault parameter.
It is **not new physics**, and any "new mechanism" claim should be dropped. What remains open
is narrower: whether an electronically tunable, moving fault in an acoustic superlattice
(items 1–4) has been described or used. The closest static precedent is Manukyan (2014). The
dynamic version was not found in this search, which is not proof of absence.

**Needed before any claim (user / library).**
- Full texts: Authier & Simon (1968), Authier's monograph (2001) chapter on planar defects,
  Manukyan (2014), Polikarpov (1987/88).
- Then compare the Entry 39 C table with the plane-wave fault formula at alpha_f = n dphi.
  This is a single analytic check with no scan.
- A targeted search for "phase-shifted" or "phase-switched" ultrasound drive in X-ray
  diffraction (two transducers or a phase-keyed drive).

**Sources (links as found).**
- Authier & Simon 1968: https://journals.iucr.org/a/issues/1968/01/00/index.html
- Kato, Usami & Katagawa 1967: https://www.cambridge.org/core/journals/advances-in-x-ray-analysis/article/abs/xray-diffraction-image-of-a-stacking-fault/03B647BC8A173C9D53A3F16B07A579E0
- Authier & Patel 1975: https://onlinelibrary.wiley.com/doi/10.1002/pssa.2210270125
- Manukyan 2014: https://www.osti.gov/etdeweb/biblio/22283837
- Polikarpov 1987: https://onlinelibrary.wiley.com/doi/abs/10.1002/pssb.2221430102
- Zolotoyabko, standing acoustic waves: https://www.osti.gov/etdeweb/biblio/20552308

---

## 2026-10-02 — Entry 56: formula-level check of the Entry 39 phase wall against Authier & Simon (1968) and Manukyan (2014)

**Sources (full text, supplied by the user; PDFs in the repo root, not committed).**
- A. Authier & D. Simon, Acta Cryst. A 24, 517–526 (1968), "Application de la théorie
  dynamique de S. Takagi au contraste d'un défaut plan en topographie par rayons X. I. Faute
  d'empilement". Takagi equations; plane- and spherical-wave cases. Sec. 4: one planar fault
  parallel to the surfaces, with the same h on both sides. chi_h^II = chi_h^I e^{i delta},
  delta = 2 pi h.u. Plane-wave amplitude eq. (4-3); its Fourier form eq. (4-5); spherical-wave
  fringes I_1, I_2, I_3; integrated intensity for a fault near a face. The case of a
  different h (moiré) is deferred to part II, which was not read.
- H. M. Manukyan, Izv. NAN Armenii Fizika 49, 360–365 (2014). Dynamical diffraction by a
  heteroepitaxial superlattice with a stacking fault at depth N_1 z_0. Satellite m is treated
  as a reflection with chi_hm = M_m chi_h. Eq. (11) gives the relative intensity change at the
  satellite centre: R_r = sin^2(alpha/2) sin(2A_1m) sin(2A_2m) / sin^2(A_m). The fault phase
  alpha = 2 pi h.u is **the same for every satellite m**.

**Mapping.** Satellite n is treated as a two-beam crystal with a = pi alpha |J_n(S)| and
thickness 1. The acoustic phase wall dphi at depth zeta_b acts as a fault with
delta = n dphi. Authier–Simon variables map as B l -> a, B eta_AS -> eta/2,
u l -> s = sqrt(a^2 + eta^2/4), l_I = zeta_b, l_II = 1 - zeta_b.

**Eq. (4-5) as used** (normalised so that delta = 0 gives the perfect crystal):

    A = a e^{i delta/2} { cos(delta/2) sin s / s - i sin(delta/2) sin(s (l_I - l_II)) / s
                          + (eta/2)/s^2 sin(delta/2) [cos(s (l_I - l_II)) - cos s] }

Derivation check. Re-derived by hand from the two-segment transfer matrix. Numerically equal
to expm products to 7e-16 (`explore/nb/e55_fault_check.py`).

Negative note. My first transcription put the last term inside the -i bracket, with a second
sin(delta/2). It matched only at delta = pi and failed elsewhere (|dI| = 0.3). It was
corrected before any comparison was recorded.

**Comparison with Entry 39 (C: wall at 0.5, ratio integral(dphi) / integral(0); window
±0.45 q, step 0.1 as in `e39_walls.py`).**

| alpha, n | a | dphi = pi/4 | pi/2 | 3pi/4 | pi | 3pi/2 | abs. integral, dphi = 0 |
|---|---|---|---|---|---|---|---|
| 0.5, 1 | 0.876 | 1.019 / 1.019 | 1.070 / 1.063 | 1.119 / 1.108 | 1.132 / 1.127 | 1.065 / 1.063 | 3.70 / 3.697 |
| 2, 1 | 3.506 | 1.069 / 1.051 | 1.198 / 1.175 | 1.287 / 1.299 | 1.327 / 1.350 | 1.168 / 1.175 | 10.1 / 10.12 |
| 2, 2 | 1.458 | 1.214 / 1.221 | 1.456 / 1.441 | 1.247 / 1.221 | 1.002 / 1.000 | 1.453 / 1.441 | 6.49 / 6.375 |
| 2, 3 | 0.383 | 1.035 / 1.011 | 1.035 / 1.007 | 1.000 / 1.002 | 1.094 / 1.013 | 1.153 / 1.007 | 0.85 / 0.873 |

Each cell is Entry 39 TT / plane-wave fault formula.

Cz (dphi = pi, n = 1), ratio pi wall / no wall:
- alpha = 2: formula 2.59 at zeta_b = 0.3 and 0.7, 1.35 at 0.5. Entry 39: 2.58 and 1.33.
- alpha = 0.5: formula maximum 1.16. Entry 39: 1.17.

Centre: |A(eta = 0)| / a = |sin(a (1 - 2 zeta_b))| holds exactly. This is the Entry 39 "echo"
law Theta = alpha_1 |1 - 2 zeta_b|, an identity of eq. (4-5) at delta = pi.

Manukyan eq. (11) follows from the same amplitude at eta = 0, using
sin^2 X - sin^2 Y = sin(X + Y) sin(X - Y):

    |A|^2 / a^2 = sin^2 a - sin^2(delta/2) sin(2 a zeta_b) sin(2 a (1 - zeta_b))

So the two papers agree with each other and with Entry 39 at the satellite centre.

**Result.**
- For n = 1 and 2, the full acoustic TT integrals of Entry 39 are reproduced by the 1968
  plane-wave stacking-fault formula with delta = n dphi, within 1–3 %.
  - This includes the absolute integrals, the odd/even selectivity, the depth dependence
    (x2.58) and the mid-depth suppression.
  - The residual 1–3 % fits what is neglected in the formula: neighbouring satellites,
    micromotion (single acoustic phase, R = 0.05), and the static-tail subtraction. Not
    separated further.
- For n = 3 (a = 0.38, nearly kinematic) the formula predicts <= 1.3 % change, while TT shows
  up to 15 %. In that small integral (0.85) this is not the fault effect. It is a different
  contribution, probably tails of neighbouring orders or the single-phase micromotion; not
  identified.
- **Verdict on Entry 39 C:** known physics, quantitatively. The acoustic phase wall is a
  stacking fault for each satellite, with fault phase n dphi. Any "new mechanism" claim for C
  is withdrawn. The B (q-wall) comparison needs Authier–Simon part II (different h, moiré
  fringes), which was not read.

**What differs from the two papers (statement of fact, not a novelty claim).**
1. In Manukyan, and in a crystallographic fault, delta = 2 pi h.u is the same for all
   satellites. The acoustic wall gives delta = n dphi, which depends on the order. The
   crystallographic counterpart is a phase slip of the modulation itself (a discommensuration
   or phase soliton of a modulated structure). That literature was not checked.
2. The fault phase and position are set by the drive: continuous in situ, and moving at the
   sound speed. The physics does not change; the experimental control does.

Whether item 1 or 2 enables a measurement that is not otherwise possible is not shown. No
calculation of such a use was made.

**Reproduction.** `python explore/nb/e55_fault_check.py`. It uses no TT solver; the Entry 39
numbers are typed in from the Entry 39 tables.
