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
