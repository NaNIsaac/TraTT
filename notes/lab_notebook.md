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

---

## 2026-10-02 — Entry 57: Chalise et al. (2026) — what "spectral weight" means, and whether our code can represent their case

**Scope (user brief).** New direction: when can the time-frequency content of a measured
diffracted intensity be read as the frequency content of the acoustic wave? Reference paper:
Chalise et al., arXiv:2603.28071v1 (30 Mar 2026), "Dynamical diffraction formalism for imaging
time-dependent diffuse scattering from coherent phonons with Dark-Field X-ray Microscopy".
Only steps 1–2 of the agreed minimal procedure were done: the definitions, and the code
capability. **No calculation was run.**

**Evidence level.** arXiv HTML text, read through a fetch tool that returns quotes and
summaries, plus the three figure images (downloaded and inspected). Numbers read off the
figures are approximate.

**1. What the paper computes.**
- Geometry: symmetric Bragg, Si (400), laser-excited Au transducer, longitudinal packet
  propagating into the sample normal to the surface (their Fig. 1a).
- Takagi–Taupin equations: Klar–Rustichelli Riccati form for X = X_1 + i X_2, R = |X(0)|^2.
  Strain enters only through the deviation parameter
  y(A) = pi V_c sin(2 theta_B) / (lambda^2 r_e f'(psi)) [Delta theta + eps(A) tan theta_B] - f'(0)/f'(psi),
  where eps(A) is the axial strain. Integration starts from the extinction depth or the back
  surface.
- Strain: from udkm1Dsim. Fig. 2a/c show bipolar decaying packets, peak |strain| ~ 1e-3
  (read off the axis), with spectral peaks at ~87 and ~113 GHz.
- Procedure: R(t) is recorded at each Delta theta and Fourier transformed at each Delta theta.
  The frequency of maximum weight vs Delta theta is linear (Fig. 1b: ~350 GHz at 53 urad to
  ~1000 GHz at 153 urad, i.e. ~6.5 GHz/urad from the plot).
- **"Weight"**: the magnitude (not the square) of the Fourier transform of R(t). Fig. 2
  compares it with the Fourier magnitude of the applied strain. **Both curves are normalised
  to their own maximum** ("Normalized Amplitude").
  - So Fig. 2 tests the *shape* of the spectrum, not its absolute scale or its linearity in
    strain amplitude.
  - The text says the weight "scales linearly with phonon population". The figure compares
    strain amplitude spectra, not populations (|amplitude|^2). The wording is inconsistent.
    Not resolved from the text.
- Stated claim: "the weight of a certain frequency omega in a phonon spectrum is exactly
  reproduced by the weight of the intensity oscillation at the same frequency".
- **Visible in their own Fig. 2.** The curves agree at the peak, but not in the tails:
  - below ~40 GHz the reflectivity spectrum lies well below the strain spectrum (~0.05 vs
    ~0.2 at ~20 GHz, both panels);
  - near 150–170 GHz there are small mismatches, including a feature in the reflectivity
    curve where the strain spectrum has a dip.
  - The paper does not discuss either. (Read off a figure; not quantified.)
- Not found in the text: the X-ray energy, the time step and window of R(t), the angular
  range used for Fig. 2, and how the angle-resolved FFTs are combined into one "frequency
  spectrum of the reflectivity".
- Also not given: a small-strain condition, an amplitude test, or any regime where the
  correspondence fails. Their resolution analysis (Sec. IV: energy bandwidth, divergence,
  objective acceptance give an oscillation damping time ~ 1/(pi FWHM)) concerns the
  instrument, not the diffraction response.

**2. Mapping to our variables (estimates, to be checked).**
- Phase-modulation index S = G u_0. For one component, u_0 = eps_0 / K, so
  S = eps_0 G / K = eps_0 G Lambda_ac / (2 pi).
  - Example: Si (400), G = 2 pi / 1.358 Å, Lambda_ac = v / nu = 8433 m/s / 87 GHz = 97 nm,
    eps_0 = 1e-3 gives **S ~ 0.7**.
  - This is a single-component estimate. For a broadband packet, S is not defined per
    component; the relevant quantity is the peak of phi = G u(z).
  - At S ~ 0.7 the first-order Bessel response is no longer linear: J_1(S) / (S/2) = 0.94.
    J_0, J_2 and the mixing between components are no longer negligible.
- Low frequencies. The ratio R = Lambda_ac / Lambda_P (Entries 42–53) corresponds here to the
  Darwin width in frequency units divided by nu.
  - Using their Fig. 1b slope (~6.5 GHz/urad) and a Darwin width of a few urad gives an
    estimate of ~10–40 GHz. **The energy is not stated, so this is an order estimate only.**
  - This puts their observed low-frequency deficit at the scale where our map leaves the
    independent-satellite regime (R >~ 0.3).
  - This is a hypothesis, not a result.
- Candidate biases for "intensity weight = strain spectrum":
  - (i) dynamical response near the Bragg peak, below a few Darwin widths in frequency;
  - (ii) strain-amplitude nonlinearity (Bessel saturation, harmonics, sum and difference
    frequencies) at S ~ O(1);
  - (iii) the normalisation hides any frequency-independent scale error by construction.

**3. Can the existing code represent their case?**

| requirement | status in this repo |
|---|---|
| symmetric Bragg, arbitrary depth profile | yes: `tratt.core.solve_bragg` (Riccati, backward from R(1) = 0), phase phi(zeta) arbitrary callable; strain enters as d phi / d zeta, phi = G u, equivalent to their y(A) eps term |
| semi-infinite crystal | approximate with a thick slab (alpha >> 1); check convergence in alpha |
| absorption | **missing** in Bragg (alpha and eta are real; the Laue gauge argument does not apply in Bragg). Needs complex eta (one dtype/validation change), or must be stated as non-absorbing |
| moving strain packet | yes: phi(zeta, t) = G u(zeta - v t) as a callable per time (`PhaseDrive` with `theta_fn` / `amp`, or a plain function) |
| udkm1Dsim strain | not available; use an analytic bipolar decaying packet with their Fig. 2 spectral peak. **Not the same profile** |
| R(t) at fixed Delta theta, FFT | yes: loop over t with `rocking(..., geometry="bragg")`, then `numpy.fft` |
| instrument resolution (Sec. IV) | not implemented; a convolution over eta in post-processing; not needed for the Fig. 2 comparison |
| absolute X-ray parameters | not in the paper text found (energy); work dimensionlessly: R (Darwin / nu) and S are the only parameters that matter for the question |

**Proposed reference case and single-variable test (not started).**
- Dimensionless, symmetric Bragg, thick non-absorbing slab, bipolar packet.
- Reference: nu_peak >> nu_Darwin and S << 1. Here intensity weight and strain spectrum
  should agree in shape, which is the paper's claim. This is the reproduction check.
- Then change one variable:
  - (a) S at fixed spectrum, to test for amplitude bias, harmonics and mixing; or
  - (b) the spectrum centre toward the Darwin scale, to test the low-frequency deficit.
- Compare un-normalised and normalised spectra.
- Stop criteria (user):
  - Stop if a bias appears only at unrealistic S or R.
  - Stop if it is fully explained by the kinematic Bessel expansion (known).
  - Stop if it only restates the paper's own limits.
  - Continue only if a usable correction or discrimination criterion (e.g. the harmonic
    ratio as an S estimator) follows at reachable conditions.

**Sources.**
- https://arxiv.org/abs/2603.28071
- https://arxiv.org/html/2603.28071v1 (figures `Untitled-1.png`, `fft.png`, `toy.png`)

---

## 2026-10-02 — Entry 58: first reference case for Chalise et al. — invalid, because the finite non-absorbing slab imprints thickness fringes

**Setup.** Script `explore/nb/e58_chalise_ref.py ref`; data in `e58_ref.json` and `e58_ref.log`.
- Dimensionless symmetric Bragg (`solve_bragg`), alpha = 2, thickness 1, non-absorbing.
- Packet phi = S exp(-tau / tau_d) sin(q0 tau), with tau = t - zeta, q0 = 2 pi 40, 3 cycles.
  The Darwin ratio is 2 pi alpha / q0 = 0.05.
- S = 0.02 and 0.05, both signs of eta. Time window T = 0.6 with 1200 steps.
- W(eta) = |FFT R(t)| taken at omega = |eta|.

**Result.**
- W peaks at eta = 242 (q0 = 251). Doubling nt and n_z changes W by 3e-4.
- The shape does not match any candidate. The maximum normalised deviation from |eps^|,
  |phi^| and |eps^|/omega^2 is the same, 0.666, for all three. This happens because the
  deviation is not set by the reference: W follows the static reflectivity R0(eta) of the slab.
- R0 has back-face thickness fringes, period ~2 pi in eta. They swing R0 by 1–3 decades:
  e.g. 3.4e-4 at eta = 220, 2.1e-3 at 256, 6.5e-5 at 327.
- The heterodyne weight scales as |r0(eta)| times the sideband amplitude, so W inherits the
  fringes. W / sqrt(R0) is smoother but still has gaps at the zeros of r0.

**Verdict.**
- This reference case is **invalid as a reproduction**.
- A finite, non-absorbing slab is not the paper's crystal. Theirs is absorbing and
  effectively semi-infinite, with integration starting from the extinction/penetration depth.
  The observed structure is an artefact of my setup, not a bias of the paper's procedure.
- One structural point does survive. The intensity oscillation at detuning eta is a
  heterodyne with the static reflection amplitude r0(eta). So W/W_max can equal the strain
  spectrum only if |r0(eta)| times the sideband transfer factor is flat over the band, or is
  divided out. In a semi-infinite crystal the tail is ~1/eta, and the transfer factor of a
  strain component (phi^ = eps^/omega) is another 1/omega. Which combination the paper's
  normalisation uses is the open definition question of Entry 57.

**Needed next (decision).** Remove the back face:
- either add absorption to the Bragg solver (complex eta; a one-line dtype change plus a test
  against the absorbing Darwin–Prins curve);
- or use a much thicker slab with an absorbing layer, which needs the same change.
Without absorption, no finite slab avoids the fringes. Not done yet.

---

## 2026-10-02 — Entry 59: absorbing reference case for Chalise et al., and a closed-form transfer function from strain spectrum to Fourier weight

**Why.** Entry 58 was invalid because a finite non-absorbing slab imprints back-face fringes.
Normal absorption was added to the Bragg solver: `solve_bragg(..., kappa)` with
eta -> eta + i kappa, tested against the Darwin–Prins root, 53 tests pass. Commit be2343d.

**Setup.**
- Scripts: `explore/nb/e58_chalise_ref.py abs darwin` (full TT) and
  `explore/nb/e59_transfer.py [fig] [weak]` (closed form, figure, weak-absorption check).
- Crystal: alpha = 8, so the Darwin half-width is 2 pi alpha = 50; kappa = 6, so
  kappa / 2 pi alpha = 0.12 and the back face is suppressed by e^-6.
- Packet as in Entry 58. Linear regime: S = 0.02 and 0.05 give the same normalised shapes.
- Darwin ratio R_D = 2 pi alpha / q0 = 0.05, 0.2, 0.5.
- The weight is now the Fourier component exactly at omega = |eta|. The nearest-FFT-bin pick of
  Entry 58 added ~15 % leakage noise.

**1. Kinematic tail (R_D = 0.05, eta > 3 Darwin half-widths).**
- |r0| eta = 25.2 = pi alpha, i.e. the tail is r0 = pi alpha / eta.
- Raw W follows |eps^| / omega^2: max normalised deviation 0.30, mostly at the low end;
  1–3 % above 0.7 q0.
- W / |r0| follows phi^ = eps^ / omega: max deviation 0.047.
- So the Fourier weight is **not** the strain spectrum. It is the strain spectrum times
  |r0(omega)| / omega ~ 1 / omega^2, because the oscillation is a heterodyne with the static
  tail amplitude.
- All spectra peak at the same frequency for this narrow packet (0.99 q0). A normalised peak
  position therefore cannot reveal the weighting; the band shape can.

**2. Closed form (first order in S).**
Linearise the Riccati equation around the static semi-infinite root R0, where
c R0^2 + i eta_c R0 + c = 0, c = i pi alpha and eta_c = eta + i kappa:

    d delta/d zeta = i c phi (1 - R0^2) - Lambda delta,   Lambda = i eta_c + 2 c R0
    int delta(0,t) e^{iWt} dt = i c S (1 - R0^2) G(W) / (Lambda + iW),   G = FT of the phase packet g
    W_pred(eta) = | conj(R0) A(|eta|) + R0 conj(A(-|eta|)) | / T,   A(W) = i c S (1 - R0^2) G(W) / (Lambda + iW)

Kinematic limit: |R0| ~ pi alpha / eta and Lambda - i eta ~ -kappa, which gives
W_kin = |R0| |G| pi alpha S / (kappa T).

| R_D | W_TT / W_pred, where W > 2 % of max | kinematic law off by up to |
|---|---|---|
| 0.05, S = 0.02 | median 0.985, range 0.971–1.010 | 0.84 |
| 0.05, S = 0.05 | median 0.984, range 0.939–1.056 | 0.90 |
| 0.2 | median 1.021, range 0.969–1.048 | 45x |
| 0.5 | median 0.993, range 0.767–1.065 | 48x |

- Time window: T = 0.6 -> 1.2 moves the median from 1.021 to 0.999, so the residual offset
  is the finite window.
- The remaining ±3.5 % scatter is not identified (candidates: the strain step at the packet
  front, the finite slab).
- **No fit parameter.** The agreement holds inside the Darwin plateau too, where W drops
  100-fold below the kinematic law.

**3. One curve.** W / W_kin, plotted against eta / (2 pi alpha) = phonon frequency / Darwin
half-width, collapses onto a single curve for all three packet spectra
(`figures/chalise_transfer.png` b).
- At kappa / 2 pi alpha = 0.12 the curve reads 0.005 (0.3 Darwin half-widths),
  0.11 (1), 0.40 (2), 0.76 (5), 0.92 (10).
- The suppression goes as ~(omega / omega_D)^2 inside the plateau.
- The factor depends on the crystal (omega / omega_D and kappa / 2 pi alpha), not on the strain
  spectrum. It is therefore a transfer function a measurer can divide out.

**4. Dependence on absorption** (closed form only, since TT was validated at 0.12; check below).

| kappa/2 pi alpha | 0.3 | 1 | 2 | 5 | 10 Darwin half-widths |
|---|---|---|---|---|---|
| 0.005 | 0.000 | 0.001 | 0.019 | 0.049 | 0.099 |
| 0.02 | 0.000 | 0.010 | 0.074 | 0.194 | 0.370 |
| 0.05 | 0.001 | 0.035 | 0.182 | 0.442 | 0.704 |
| 0.12 | 0.005 | 0.109 | 0.396 | 0.756 | 0.919 |
| 0.3 | 0.021 | 0.280 | 0.679 | 0.929 | 0.981 |

- Weak absorption extends the suppression far into the tail. There the probed depth is set by
  the tail extinction depth ~ eta / (pi alpha)^2 rather than by absorption: in
  Lambda - i eta = -kappa + 2 c R0, the term |2 c R0| ~ 2 (pi alpha)^2 / eta competes with kappa.
- Order estimate (not computed for a real crystal): for Si reflections at ~10 keV,
  kappa / 2 pi alpha is of order 1e-3–1e-2. That would make the correction large across the
  band of their Fig. 2.
- TT check at 0.024 (alpha = 40, kappa = 6, R_D = 0.2; `e59_transfer.py weak`): W_TT / W_pred
  median 0.993, range 0.944–1.044. The closed form holds at weak absorption too.

**5. Relation to Chalise et al. (Entry 57).**
- Their Fig. 2 normalises both spectra to their peaks and reports exact overlap at the peak.
- In our model a peak-normalised comparison of a narrow band also overlaps at the peak.
  Off the peak the weight carries |r0| |phi^| times the transfer factor.
- Their visible low-frequency deficit is consistent in sign with the near-Darwin suppression of
  section 3.
- Their weighting of the strain curve, energy and absorption are not given, so no quantitative
  comparison was made. Their model is not claimed to be wrong.

**Verdict.**
- The Fourier weight of R(t) equals the strain spectrum only after dividing out a transfer
  function H(omega) = W / |eps^|.
- H(omega) is given in closed form at first order in S by the static root R0(omega). Its
  ingredients are the X-ray reflection, the extinction and the absorption.
- It is verified against full TT to 2–5 % for kappa / 2 pi alpha = 0.12 (R_D = 0.05–0.5) and
  0.024 (R_D = 0.2).
- This is a usable correction ("divide the measured weight by H") at reachable conditions,
  which meets the continue criterion of the user's plan.
- Novelty is not assessed. Linear-response treatments of time-resolved X-ray diffraction may
  contain the same form; that literature was not searched.

**Not done.**
- Anomalous absorption (chi_h''), asymmetric geometry, instrument resolution.
- Nonlinear S (harmonics, mixing).
- A real-crystal parameter set.
- A literature check of the linear-response form.

---

## 2026-10-02 — Entry 60: literature status of the Entry 59 transfer function (search only)

**Evidence level.** Web search results and one abstract. No full text was read.

**Found.**
- Shayduk, Herzog, Bojahr, Schick, Gaal, Leitenberger, Navirian, Sander, Goldshteyn, Vrejoiu,
  Bargheer, Phys. Rev. B 87, 184301 (2013), arXiv:1303.5643: "Direct time domain sampling of
  sub-THz coherent acoustic phonon spectra in SrTiO3 using ultrafast X-ray diffraction".
  - They determine a generated phonon spectrum from ultrafast XRD and find a k-dependent
    phonon lifetime (Akhiezer 1/omega^2).
  - Quote from the abstract: "the observed shift of the synthesized spectrum to the higher q is
    discussed in the framework of non-linear effects".
  - **Hypothesis to check, not a claim:** the dynamical suppression of the weight at low
    frequency (Entry 59, section 3) also shifts an apparent spectrum to higher q. Whether their
    analysis divides out a transfer function is not known from the abstract.
- "Reconstructing longitudinal strain pulses using time-resolved x-ray diffraction"
  (ResearchGate entry; fetch refused): strain reconstruction from UXRD is an established topic.
- "Concepts and use cases for picosecond ultrasonics with x-rays", arXiv:2302.14116 (review,
  Bargheer group). Likely to state the standard kinematic relation. Not read.
- "Time-Dependent Dynamical Diffraction Theory for Phonon-Type Distortions", Springer chapter
  (doi via link.springer.com/chapter/10.1007/978-1-4615-0387-3_10). Not read.

**Status.**
- The kinematic heterodyne weighting (W ~ |r0| |phi^|, i.e. ~|eps^|/omega^2) is probably known
  in this community.
- Whether the full first-order dynamical closed form (Entry 59 section 2, with 1 - R0^2 and
  the 2 c R0 extinction term) and its one-curve suppression near the Darwin width are published
  is **open**.

**Needed.** Read Shayduk 2013 (analysis section), the 2023 review, and the Springer chapter.
Then decide whether Entry 59 is a correction worth a short paper.

**Update (same day, full text of Shayduk et al. 2013, arXiv pp. 1–6 read).**
- **Their estimator is not the Fourier weight.** They use the time-averaged intensity excess
  <I_p(Q)>_t - I_up(Q) = alpha A(q)^2 (their eq. 2). This is second order in the strain, the 1D
  analogue of thermal diffuse scattering (their eqs. 3–8).
- They name the first-order oscillation explicitly ("interference of X-rays scattered from a
  moving grating and the static component ... an X-ray analog of Brillouin oscillations") and
  average it away. Chalise et al. use exactly that oscillation, and so does our Entry 59.
- **Their validity conditions** (p. 4) include q >> 1/l, with l the X-ray extinction (or coherence)
  length; time averaging over many periods; displacements << interatomic distance; q << G.
  Their measured peaks (q ~ 0.045 and 0.09 rad/nm, STO 002 at 8 keV) sit well off the Darwin
  width. They state that the X-ray extinction is negligible there and only absorption matters.
- **Hypothesis of this entry revised: unlikely.** Their peak drift to higher q grows with time as
  the packet goes deeper (Fig. 4b). A static transfer function cannot produce a time-dependent
  drift, and their q are far from the Darwin width.
- What remains: the *first-order* (oscillation) estimator near and around the Darwin width, and
  in weakly absorbing crystals far into the tail (Entry 59 section 4). This is the regime that
  Shayduk's condition q >> 1/l excludes and that Chalise's analysis uses.
- Also relevant: Shayduk eq. (14) is the absorption decay time tau_abs = 2 mu c_s / sin theta.
  That is the same physics as kappa in our closed form.

**Update 2 — likely closest prior art for the closed form (abstract level).**
- P. Sondhauss & J. S. Wark, "Extension of the time-dependent dynamical diffraction theory to
  'optical phonon'-type distortions: application to diffraction from coherent acoustic and
  optical phonons", Acta Cryst. A 59 (2003). IUCr page refused (403); index:
  https://journals.iucr.org/a/issues/2003/02/00/index.html
- Book chapter: "Time-Dependent Dynamical Diffraction Theory for Phonon-Type Distortions", in
  *Nonlinear Optics, Quantum Optics, and Ultrafast Phenomena with X-Rays* (Springer, 2003).
- Abstract: a generalised, time-dependent Takagi–Taupin equation; "numerical and perturbative
  analytical solutions ... in the presence of such phonons", comparing acoustic and optical
  phonons.
- **Consequence.** A first-order perturbative solution for acoustic phonons in Bragg geometry
  may already exist there. The Entry 59 closed form must not be called new until this paper is
  read. What may still be new is the *use*: the weight-to-spectrum transfer function H(omega),
  its one-curve collapse in omega/omega_D, the strong dependence on kappa/2 pi alpha, and the
  comparison with the peak-normalised spectra of Chalise et al. That too is open until the
  paper is read.
- **Needed (library):** Sondhauss & Wark 2003, both versions.

## 2026-10-07 — Entry 61: how far the Entry 59 transfer function is prior art (secondary sources only)

**Why.** Entry 60 flagged Sondhauss & Wark, Acta Cryst. A59 (2003) as the likely closest prior
art for the Entry 59 closed form and said the paper must be read before any novelty claim. The
paper is still closed access. This entry does not read it; it establishes what the surrounding
literature does, which bounds the claim from the other side.

**Evidence level.** One reconstructed abstract, two full-text PDFs read by keyword, two
bibliographic APIs. The primary source is still unread.

**1. The paper itself.**
- OpenAlex (`api.openalex.org/works/doi:10.1107/S010876730201704X`): 43 citations, no open-access
  location, abstract reconstructed. Semantic Scholar reports 1 citation for the same DOI and no
  abstract, so its count is an indexing gap, not a measure of impact. Use OpenAlex here.
- Abstract confirms what Entry 60 assumed: "Numerical and perturbative analytical solutions of
  generalized Takagi-Taupin equation are presented in presence of such phonons." So a first-order
  perturbative solution for phonon distortions very likely exists in that paper.
- Access attempts that failed: IUCr (403), PubMed (cookie wall), Springer chapter (paywall),
  Semantic Scholar (abstract elided by publisher), Elsevier for the companion paper
  Sheppard, Sondhauss, Merlin, Bucksbaum, Lee & Wark, Solid State Commun. 136, 181 (2005),
  "Simulations of the phonon Bragg switch in GaAs", which uses the same perturbative approach.
  The Lund repository record for that one carries no full text either.

**2. What the picosecond-ultrasonics community does (arXiv:2302.14116, Bargheer group review,
29 pp, full text read by keyword).**

| term | occurrences |
|---|---|
| Takagi | 0 |
| Darwin | 0 |
| extinction | 0 |
| transfer function | 0 |
| Sondhauss | 2 |

- Both "Sondhauss" hits are unrelated co-authorships in the reference list (InSb acoustic phonons,
  Appl. Phys. A 75, 467 (2002); a streak camera, Opt. Lett. 35, 3219 (2010)). **The review does not
  cite Sondhauss & Wark 2003.**
- "dynamical" appears three times, twice for the numerical `udkm1Dsim` libraries and once for
  peak broadening in inhomogeneously strained layers ("require explicitly modelling the Bragg peak
  evolution via kinematical or dynamical x-ray diffraction").
- So the standard route is kinematic, or numerical dynamical simulation of the whole peak. No
  analytic frequency-dependent correction is in use.

**3. What the direct comparison case does (Chalise et al. 2026, arXiv:2603.28071, full text read
by keyword).**

| term | occurrences |
|---|---|
| Takagi | 12 |
| dynamical diffraction | 9 |
| kinematic | 9 |
| extinction | 2 (both as the depth at which the iteration starts) |
| Darwin | 0 |
| transfer function | 0 |
| Sondhauss / Wark | 0 |

- They solve the Takagi-Taupin system numerically through the depth and take `R = |X(0)|^2`.
  That is the same physics as our solver, so the suppression of Entry 59 section 3 is present in
  their forward model.
- They do not write a first-order closed form, do not mention the Darwin width, and do not divide
  anything out when they report spectral weight.

**4. Citing works of Sondhauss & Wark 2003** (OpenAlex, 50 fetched). Open-access ones are
"Coherent phonon control" (2007), "Tuning and Tracking of Coherent Shear Waves in Molecular Films"
(2018), an InSb heat-conduction paper (2008) and a grazing-incidence TDS paper (2003). None
reproduces or applies a first-order Bragg transfer function. Several entries are mis-attributed
(uranium solution chemistry), so the list is noisy.

**Conclusion (bounded).**
- The first-order perturbative *solution* of the time-dependent Takagi-Taupin equation for phonon
  distortions is probably in Sondhauss & Wark 2003. Do not claim that step as new.
- The *use* of it as a transfer function H(omega) between the strain spectrum and the measured
  Fourier weight -- its collapse in omega/omega_D, its dependence on kappa / 2 pi alpha, and the
  instruction to divide it out -- has no visible precedent: the field's own 2023 review does not
  cite the paper and uses no such correction, and the 2026 paper that needs it does not have it.
- This bounds, but does not close, the question. Reading Sondhauss & Wark remains required before
  submitting the correction as a result.

**Poster.** No change. The transfer function appears only in the Outlook line, where it is stated
as a direction, not as a priority claim.

**Needed (library).** Sondhauss & Wark, Acta Cryst. A59, 7-13 (2003), DOI
10.1107/S010876730201704X -- via the Sogang subscription or an IUCr reprint request. Second
choice: Sheppard et al., Solid State Commun. 136, 181 (2005), DOI 10.1016/j.ssc.2005.05.056.

---

## 2026-10-07 — Entry 62: the quasi-static number for the poster's experimental window

**Why.** The poster called the model "time-dependent two-beam Takagi-Taupin". The solver is not
that: `core.py` integrates one frozen phase profile per instant, and the header of this notebook
says "stroboscopic". The label also collides with Sondhauss & Wark 2003, whose title is literally
"time-dependent dynamical diffraction theory" (Entry 60, 61) -- a reader from that community would
take the two to be the same object. The fix is to name the approximation and give its number.

**Parameter.** The lattice is static during one X-ray traversal when the acoustic phase advanced in
that time is small:

    beta = Omega L / c = 2 pi f L / c

with L the slab thickness (the longest path in the problem) and c the vacuum speed of light. The
stricter-looking version with the Pendellosung length, 2 pi f Lambda_P / c, is smaller, so beta is
the conservative choice.

**Numbers**, for the recommended window of Entry 20 (Si 440, symmetric Laue, L = 3 Lambda_P):

| E keV | Lambda_P um | L um | f GHz | beta = 2 pi f L / c | 2 pi f Lambda_P / c |
|---|---|---|---|---|---|
| 17.48 | 54.34 | 163 | 1.1 | 3.8e-3 | 1.3e-3 |
| 17.48 | 54.34 | 163 | 2.0 | 6.8e-3 | 2.3e-3 |
| 30 | 98.76 | 296 | 0.6 | 3.7e-3 | 1.2e-3 |
| 30 | 98.76 | 296 | 2.0 | 1.2e-2 | 4.1e-3 |

- Worst case in the window: **beta = 1.2e-2**.
- Stated as a time rather than a phase: at the thickest, fastest corner the X-ray crosses the
  crystal in L/c = 0.99 ps, which is **0.2 % of the 500 ps acoustic period**.

**Reading.** The frozen-profile treatment is not a convenience; in this window the lattice moves by
0.2 % of a period while the field is inside. The poster should say so and drop "time-dependent".

**Not covered by this number.** The X-ray response time of the crystal itself (the field needs
~Lambda_P / c to settle into the dynamical state, the regime of Shvyd'ko & Lindberg) is a separate
question from the acoustic motion, and the stroboscopic model does not address it. Entry 21's
caveat stands: any genuine dependence on v_s beyond t_tr is outside these maps by construction.

---

## 2026-10-08 — Entry 63: what a surface acoustic wave has to do with this

**Why.** SAW X-ray diffraction is the largest experimental community in acoustically driven
diffraction, and the poster says nothing about it. Entry 21 only recorded that Zolotoyabko's
standing-SAW work is a periodic steady state. The geometric relation was never worked out, so
"why not just do it with a SAW?" has had no answer. Script: `explore/nb/e63_saw.py`.

**1. The geometries are different, and the solver cannot bridge them.**

| | this work | SAW (Rayleigh) |
|---|---|---|
| acoustic wavevector | along the depth | in the surface |
| amplitude behind the front | uniform | decays into the depth |
| what grows with time | the modulated **depth** | the modulated **lateral extent** |
| satellite displacement | in the rocking variable | needs the in-plane component (hence "sagittal") |

The solver is one-dimensional in depth. An in-plane acoustic wavevector changes the in-plane
momentum transfer and cannot be represented in it at all. **No result of this work transfers to
SAW geometry as it stands.** That is a statement about the model, not about the physics.

**2. What does carry over: the law is about the depth profile, not about the front.**

The accumulated-coupling law (Entries 5, 28, 30) reads the depth profile through a Bessel kernel,
`Theta_n = int_0^1 J_n(S A(zeta)) d zeta`, and Laue displays it as `sin^2(pi alpha |Theta_n|)`.
Nothing in that derivation requires `A` to be a travelling step. Tested here with a **static
exponential skin**, `A(zeta) = exp(-zeta / zeta_d)`, carrier still travelling, averaged over 8
acoustic phases (alpha = 3, alpha_ac = 20, n = 1, n_z = 4096):

| zeta_d | S = 0.5 | S = 1.5 | S = 3.0 |
|---|---|---|---|
| 0.02 | 8.6e-3 | 5.0e-3 | 2.9e-3 |
| 0.05 | 8.3e-3 | 3.9e-3 | 2.6e-3 |
| 0.10 | 6.7e-3 | 4.2e-3 | 1.1e-2 |
| 0.30 | -3.4e-2 | -7.2e-2 | 4.0e-2 |
| 1.00 | -1.6e-1 | 2.8e-2 | -2.9e-3 |

(absolute difference, full TT minus law). **For skins thinner than 0.1 L — which is the whole SAW
regime, see section 3 — the law holds to <= 1.1e-2 absolute.** It degrades for a thick skin, where
the drive is strong over the whole depth and the known `(alpha/alpha_ac)^2` residual plus the
phase average take over. This is a genuine extension of the law's domain: step onset, finite ramp,
and now a static depth profile.

**3. How much of the extinction depth a SAW can modulate — a bound, not an estimate.**

A Rayleigh wave is two partial waves decaying as `exp(-kappa_{L,T} z)` with
`kappa_{L,T} = q sqrt(1 - c_R^2 / c_{L,T}^2)`, `q = 2 pi / Lambda_ac`. The slower one sets the
skin. For Si (`c_L = 8435`, `c_T = 5846` m/s from c11, c44):

| c_R | kappa_T / q | skin |
|---|---|---|
| 4920 m/s | 0.540 | 0.295 Lambda_ac |
| 5378 m/s (0.92 c_T) | 0.392 | 0.406 Lambda_ac |

So the fraction of the extinction depth that a SAW modulates is

    skin / Lambda_P = (0.30 to 0.41) * R,        R = Lambda_ac / Lambda_P

— **frequency- and material-independent in R**. The resolved-channel condition of Entry 54 is
`R <= 0.34-0.57`, so

    modulated fraction <= 0.10 to 0.17.

A SAW cannot fill the extinction depth and keep its satellites resolved. The two requirements pull
in opposite directions: resolving channels wants a short acoustic wavelength, filling the depth
wants a long one. Si numbers (`c_R = 4920`, Lambda_P from Entry 20):

| f GHz | Lambda_ac um | skin um | Si 440 17 keV: R / frac | Si 440 30 keV: R / frac |
|---|---|---|---|---|
| 0.2 | 24.6 | 7.25 | 0.453 / 0.133 | 0.249 / 0.073 |
| 0.5 | 9.84 | 2.90 | 0.181 / 0.053 | 0.100 / 0.029 |
| 1.0 | 4.92 | 1.45 | 0.091 / 0.027 | 0.050 / 0.015 |
| 2.0 | 2.46 | 0.73 | 0.045 / 0.013 | 0.025 / 0.007 |
| 5.0 | 0.98 | 0.29 | 0.018 / 0.005 | 0.010 / 0.003 |

**4. What a skin reads as a filled depth.** With `t_eff = Theta_n(skin) / J_n(S)`, a skin of depth
`zeta_d` reads like a filled depth of

| S | t_eff / zeta_d (n = 1) |
|---|---|
| 0.5 | 1.02 |
| 1.5 | 1.23 |
| 3.0 | 3.09 |

The ratio exceeds 1 at large S because the decaying tail of the skin sweeps down through the
maximum of `J_1` (at S ~ 1.84): depth where the amplitude has fallen contributes *more* coupling
than the surface. A deep crystal can therefore still show the oscillation — but only by changing
frequency or thickness, which is a steady-state measurement.

**5. The decisive point: a SAW has no depth-filling transient.**

The skin depth is set by the acoustic frequency and does not grow. Once the surface wavefront has
passed a point, the depth profile under it is established and static. The transient in this work
exists because the acoustic front travels **into** the depth that the X-ray samples; for a SAW the
front travels **across** the surface, and the depth profile is slaved to it.

A SAW does have a transient, but it is a different one: the lateral filling of the X-ray
footprint, `footprint / v_R` — about 20 ns for a 100 um footprint, coincidentally the same order
as our `t_tr`. That produces an inhomogeneously illuminated surface, not a depth-graded crystal,
and it is not what the filled-depth law describes.

**Consequence for the poster.** The experiment this work predicts is a **bulk** wave launched into
the depth from a transducer, and shear-polarised so that `h . u != 0` in symmetric Laue (Entry 6,
Entry 20, handoff 7.4). That is already the poster's geometry; it should be said out loud, because
a visitor from the SAW community will otherwise assume a surface wave.

**Not done.**
- No calculation with an in-plane acoustic wavevector; the solver cannot do one.
- The exponential skin is a model. The real Rayleigh profile is a sum of two partial waves whose
  normal displacement changes sign with depth, and the `h . u` projection varies with depth too.
- The 2017 sagittal-SAW Bragg paper (handoff ref [3]) has still not been read.
- Nothing here says a SAW experiment is uninteresting; it says it measures a different quantity.

---

## 2026-10-14 — Entry 64: Sondhauss & Wark (2003) read — what it already contains, and what Entry 59 adds

**Source.** P. Sondhauss & J. S. Wark, Acta Cryst. A59, 7-13 (2003), DOI 10.1107/S010876730201704X.
Full text, all 7 pages, read. File: `SondhaussWark2003_ActaA59_7.pdf` in the repo root, untracked
(copyright), next to the Authier & Simon PDF. The page range cited in Entries 60-61 and the handoff
(100-107) was wrong and has been corrected to 7-13.

**1. What the paper does.**
- Sec. 2-3: generalised time-dependent TT equation (their eq. 14) for non-uniform (optical-phonon)
  displacements; the classical time-dependent TT equation is recovered as the special case of uniform
  displacement (eq. 16). The explicit time-derivative terms matter only when Omega/omega is not
  negligible (polaritons); for acoustic phonons they shift satellites by eq. (36), a relativistic
  Doppler effect. Consistent with Entry 62 (quasi-static for our window).
- Sec. 4, the relevant part: **first-order perturbative solution.** Split chi'_h = chi_h + delta chi_h
  cos(Omega t - q.r) and D'_h = D_h^(p) + delta D'_h, drop delta chi * delta D. The result (eq. 25) is a
  linear system with constant coefficients and the phonon as the source term. Two-beam case (eq. 26),
  with the single-mode plane-wave ansatz delta D'_h = d_h^+ e^{i(Omega t - q.r)} + d_h^- e^{-i(...)}
  (eq. 28), gives an algebraic solution (eq. 29):

      d_h^± = (C/2) [chi_h dchi_hbar C D_0^(p) ± A_0 dchi_h D_h^(p)] / [±A_0 (±A_h + 2 beta_h) - chi_h chi_hbar C^2]

  with A_0, A_h = (2/k) s.q - 2 Omega/omega (eq. 30) and beta_h the deviation parameter (eq. 27, 33).
  Total field = D_h^(p) + satellites at |G_h ± q| oscillating at Omega (eq. 31).
- Verified against full numerics for InSb 111 Bragg, 10 um, |q| = 1e8 rad/m, strain 0.1 % (Fig. 3:
  "difference ... cannot be seen on the scale of the graph").
- **Stated limits** (p. 12):
  - reasonable while satellite peak intensity <= 10 % of the main peak;
  - **not safe when the sidebands are within about one rocking-curve width of the main peak**, i.e.
    phonon wavelength comparable to or longer than the extinction depth; Fig. 4 (|q| = 5e6 rad/m,
    strain 0.01 %) shows visible disagreement there;
  - inaccurate for crystals thin enough to show Laue (thickness) oscillations.
- Proposed use (p. 12, one paragraph, no demonstration): because eq. (25) is linear, mode
  contributions superpose, so a least-squares fit of a sum of analytic solutions to a measured
  rocking curve "directly reveals information about the modes' occupation numbers", with an
  ambiguity from the missing phase.

**2. Against Entry 59, point by point.**

| question (from Entry 60-61) | Sondhauss & Wark 2003 | Entry 59 |
|---|---|---|
| first-order perturbative solution of TT for an acoustic phonon in Bragg | **yes** (eq. 25-29) | yes, independently, in Riccati form |
| form | single mode, plane-wave ansatz with constant d^±, two-field algebraic system | arbitrary packet G(W), linearised Riccati about the semi-infinite root R0 with boundary condition at depth, Lambda = i eta_c + 2 c R0 |
| absorption | through complex chi_0 in beta_h (implicit) | explicit kappa; the weight depends strongly on kappa / 2 pi alpha (Entry 59 sec. 4) |
| Fourier weight of the intensity oscillation as a spectral estimator | **no**; the oscillating cross term D^(p)* d^± is present in eq. 31 but not discussed; Fig. 3-4 show rocking curves | yes: W = heterodyne with r0, the transfer function H(omega) |
| dividing out a crystal-dependent factor | no | yes |
| one-curve collapse in omega / omega_D | no | yes (Entry 59 sec. 3) |
| regime near the Darwin width (sidebands within ~1 rocking width) | **explicitly declared unsafe** for their solution (Fig. 4) | verified against full TT at R_D = 0.2, 0.5, i.e. inside and at the edge of the plateau: W_TT / W_pred median 0.99-1.02 |
| use for measuring phonon populations | proposed: fit rocking curve with a sum of mode solutions | proposed: divide the measured Fourier weight by H(omega) |

**3. Verdict.**
- The first-order perturbative solution itself is **prior art**. Entry 59 must cite it and must not
  call the linearisation new.
- The idea of using the linear solution to extract phonon mode content is also **present in kind**
  (S&W p. 12), as a rocking-curve fit, without demonstration.
- Not in S&W:
  - the time-Fourier weight of R(t) as the estimator, and its heterodyne weighting by r0;
  - the explicit transfer function H(omega), its collapse on omega / omega_D, and its kappa dependence;
  - a first-order form that stays valid where S&W say theirs fails: sidebands within about one
    rocking-curve width.
- That last point is the sharpest difference, and it is **stated, not shown**. S&W's failure in
  Fig. 4 is for a different crystal, and their eq. 29 has not been evaluated in our variables. The
  likely reason for the difference is that the constant-d ansatz ignores the depth dependence of
  D^(p) and the surface boundary condition of delta D, which Entry 59 keeps through R0 and Lambda.
  This is a hypothesis.
- Equivalence of the two first-order forms far from the Darwin width (where both should be the
  kinematic heterodyne) is plausible but not checked.

**4. Smallest test that would turn the stated difference into a shown one (run; see section 5).**
Evaluate S&W eq. 29 in the dimensionless symmetric-Bragg variables of `e59_transfer.py` (single mode,
same alpha, kappa). Compare it with full TT and with the Entry 59 form at R_D = 0.05 (both should
agree) and R_D = 0.5 (S&W expected to fail). Stop if S&W eq. 29 also holds at R_D = 0.5: then the
remaining difference is only the estimator (Fourier weight vs rocking curve) and the H(omega)
presentation.

**Other notes.**
- Their 10 % satellite-to-main-peak rule is a usable amplitude criterion for the nonlinear-S item of
  Entry 59 "Not done".
- Not read: Sheppard et al. 2005 (same approach, phonon Bragg switch); Chukhovskii & Förster 1995,
  Wark & Lee 1999 (time-dependent TT).

**5. Result of the section 4 test (run 2026-10-14).** Script `explore/nb/e64_sw_compare.py`, log
`e64_sw.log`, data `e64_sw.json` (9 s).
- Setup as in Entry 59: alpha = 8, kappa = 6 (kappa / 2 pi alpha = 0.12), S = 0.02. A single mode
  phi = S cos(W (t - zeta)) fills the crystal. The estimator is the Fourier component F of R(t) at W,
  over a rocking scan.
- TT convergence (nt x2, n_z x2) is 1.4e-5 of max |F|.
- Three first-order predictions, from the two-field form of the `solve_bragg` equation:
  - **SW** is S&W eq. 28-29 read literally: constant perfect-crystal amplitudes (1, R0), constant
    d^±, satellite amplitude read as the Dh perturbation. The 2 Omega/omega terms are dropped.
  - **SWn** is the same, normalised to the incident wave (dR = B - R0 A).
  - **Riccati** is the Entry 59 form. With the true Bloch wave (1, R0) e^{c R0 zeta} instead of
    constants, the normalised two-field solution equals it identically (max difference 6e-16 of
    max |F|). So the only modelling difference is that ansatz.

| W in Darwin half-widths (R_D) | SW at satellite eta = ±W | SWn at satellite | Riccati at satellite | SW worst near main peak, abs(eta) < 1.5 w_D | Riccati complex error, max |
|---|---|---|---|---|---|
| 20 (0.05) | 1.022 | 1.021 | 1.005 | 3.5 | 0.01 |
| 5 (0.2) | 1.220 | 1.207 | 1.003 | 3.7 | 0.01 |
| 2 (0.5) | 1.687 | 1.567 | 1.001 | 4.4 | 0.01 |
| 1 (1.0) | 3.014 | 1.673 | 1.000 | 6.2 | 0.00 |
| 0.5 (2.0) | 8.450 | 12.96 | 1.000 | 9.4 | 0.00 |

(The ratios in the satellite columns are |F_pred| / |F_TT|. The other errors are |F_pred - F_TT| / |F_TT|
where |F_TT| > 2 % of max. Away from the main peak and the satellites, SW stays within 0.04-0.19.)

**Reading.**
- S&W's own warning (p. 12, Fig. 4) is reproduced in our variables. Their solution is fine at the
  satellite when the satellite sits 20 Darwin half-widths out (2 %). It is 22 % off at 5, 69 % at 2,
  a factor 3 at 1 and a factor 8 at 0.5.
- It also fails near the main peak (factor 3-9) at every W. That is where the perfect-crystal field
  is evanescent and a constant D^(p) is worst.
- Normalising to the incident wave does not repair it. **The cause is the constant-amplitude ansatz:
  it ignores the extinction of D^(p) with depth.** With the Bloch wave restored, the solution
  becomes the Entry 59 form and agrees with TT to 1 % everywhere, including inside the Darwin plateau.
- The section 3 difference is now **shown**, at one absorption (kappa / 2 pi alpha = 0.12) and one
  small S. In the regime S&W declare unsafe, which is the regime Entry 59 sec. 3 identifies as the
  biased part of the Chalise estimator, the first-order theory still works. It must use the
  Bloch-wave depth dependence. Entry 59's Riccati linearisation does this automatically.

**Novelty status (bounded, updated).**
- Prior art: the linearisation itself (S&W 2003), and fitting mode contents with a linear solution
  (S&W, proposed).
- Not in S&W, and now demonstrated:
  - a first-order form valid within about one rocking width of the main peak;
  - its use as the transfer function H(omega) for the time-Fourier estimator;
  - the omega / omega_D collapse and the kappa dependence.
- Fair caveats:
  - this is our reading of their eq. 28-29, and their own numerics (Fig. 3-4) may have used more;
  - correcting a single-Bloch-wave ansatz is a modest step, so a referee may call it routine. The
    paper case must rest on the measurement correction, not on the linearisation.
- Still unread: Sheppard et al. 2005, Chukhovskii & Förster 1995, Wark & Lee 1999.

---

## 2026-10-14 — Entry 65: the transfer function for a real crystal (Si 400), and where it stops being linear

**Why.** Entries 59 and 64 give the crystal factor H in dimensionless form (frequency over the
Darwin half-width x, absorption g = κ/2πα). This entry puts in Chalise et al.'s crystal and asks
two things: how large is the bias across their spectral band, and is the closed form still right
there? The regime is far out in the weakly absorbing tail (x ~ 100), where nothing had been
checked against TT. Scripts: `explore/nb/e65_si400.py` (closed form, full-TT check, figure) and
`explore/nb/e65_conv.py` (depth-step convergence and the linear-regime test). Logs:
`e65_si400.log`, `e65_conv.log`. Figure: `figures/si400_transfer.png`.

**1. Mapping (sigma polarisation, X0h constants of Entry 20, symmetric Bragg).**
- Λ = λ sinθ_B/|χ_rh| (the code's Λ with γ = sinθ_B).
- g = κ/2πα = χ_i0/|χ_rh|. With μ = 2πχ_i0/λ and κ = μL/sinθ_B, the factor sinθ_B cancels.
- x = f/f_D, with f_D = v_L/Λ and v_L = 8433 m/s.

| E keV | Λ µm | f_D GHz | g |
|---|---|---|---|
| 8 | 11.33 | 0.744 | 0.0458 |
| 17.48 | 11.54 | 0.730 | 0.0100 |

So Chalise's 87–113 GHz peak sits at x ≈ 120–150, and their 20–170 GHz band at x ≈ 27–230.

**2. Crystal factor across their band (closed form, single mode at the satellite).**

| f GHz | 5 | 10 | 20 | 40 | 60 | 100 | 170 |
|---|---|---|---|---|---|---|---|
| H, 8 keV | 0.52 | 0.78 | 0.93 | 0.98 | 0.99 | 1.00 | 1.00 |
| H, 17.48 keV | 0.14 | 0.26 | 0.48 | 0.74 | 0.85 | 0.94 | 0.98 |

- **H(20 GHz)/H(100 GHz) = 0.51 at 17.48 keV and 0.93 at 8 keV.**
- A peak-normalised Fourier-weight spectrum taken at 17 keV therefore under-reads 20 GHz by a
  factor 2 after any kinematic correction. Their Fig. 2 shows a low-frequency deficit of the same
  sign; its size is about 4× (Entry 57, read off the figure).
- Their energy and their weighting of the strain curve are not stated (Entry 57), so this stays
  consistency in sign and order of magnitude, not a quantitative match.
- The factor comes from weak absorption. Where κ is small, the tail extinction 2cR0 ~ (πα)²/η
  competes with absorption out to tens of Darwin half-widths (Entry 59 sec. 4).

**3. Full-TT check: two separate effects.**
- **(a) Numerics.** At 24 RK4 steps per cycle of η, the phase error accumulates over ~10⁴
  cycles. At g = 0.01, x = 30 the ratio TT/closed goes 0.917 → 0.954 → 0.956 for 24/48/96 steps
  per cycle; at x = 100 it goes 0.680 → 0.744 → 0.746. The h⁴ scaling holds, and 48 steps per
  cycle is converged to < 0.5 %. **The `tt` run in `e65_si400.log` used 24 steps per cycle and
  is not converged at x ≥ 30. Do not quote it.**
- **(b) Physics: the satellite stops being linear.** The converged residual at S = 0.02 is not
  numerical. At S = 0.002, with 48 steps per cycle, the closed form holds:

| g | x | S | ν = linear satellite amplitude | 1 − TT/closed |
|---|---|---|---|---|
| 0.0100 | 30 | 0.02 | 0.257 | 0.044 |
| 0.0100 | 100 | 0.02 | 0.447 | 0.254 |
| 0.0100 | 30 | 0.002 | 0.026 | 0.002 |
| 0.0100 | 100 | 0.002 | 0.045 | 0.013 |
| 0.0458 | 30 | 0.002 | — | 0.001 |
| 0.0458 | 100 | 0.002 | — | 0.005 |
| 0.12 (Entry 64) | 20 | 0.02 | 0.041 | −0.005 |

- The deviation follows ν², where ν = |c S (1 − R0²) / (2(Λ − iW))| is the satellite amplitude
  the linear theory predicts.
- In words: when absorption is weak and the strain fills the penetration depth, the satellite
  builds up its own reflectivity. It reaches the regime where its coupling πα|J_1(S)| over the
  penetration depth is no longer small, and the first-order (kinematic-in-S) weight over-reads.
- This is the same physics as the rest of this work: the satellite becomes a dynamical reflection
  with its own extinction length Λ/|J_n(S)|.

**4. Validity rule for the Fourier-weight estimator, and what it means for Chalise.**
- The linear transfer function H applies while ν ≲ 0.2 (error ≲ 4 %).
- For a packet rather than a filling mode, ν scales with the strained depth that the X-ray
  samples coherently: ν ~ π|J_1(S)| ℓ/Λ with ℓ = min(packet length, penetration depth).
- Chalise's packets are ~3 cycles of ~97 nm (ℓ ≈ 0.3 µm), with |strain| ~ 1e-3, so S ≈ 0.7
  (Entry 57) and J_1 ≈ 0.33. That gives **ν ≈ 0.03, linear: H applies to their case.**
  This is an order estimate, not a TT run of their packet.
- For long acoustic trains — the GHz transducer drive of this project — ν is large by design, and
  the satellite forms dynamically (Entries 1, 9). The two halves of the poster meet here:
  - **short packet** → linear, H-corrected spectroscopy;
  - **long train** → satellite Pendellösung.
  - The crossover is the satellite's own extinction length reaching the strained depth.

**Verdict.**
- For Si 400 the crystal factor is a factor-2 effect across a realistic DFXM band at 17 keV, and
  < 10 % at 8 keV.
- The closed form is verified against converged TT to ≤ 0.5 % in the linear regime at both
  energies, out to x = 100.
- The estimator itself has a linearity limit, ν ≲ 0.2. It is a property of the measurement, not of
  the solver, and it ties the Outlook back to the formation physics.

**Not done.**
- A TT run of an actual few-cycle packet at Si 400 parameters, which would replace the ν
  estimate.
- Anomalous absorption χ_ih in Bragg.
- Their exact energy and normalisation (not stated in the preprint).
- (nothing pending)

**Update (same night).** At 8 keV, S = 0.02, converged with 96 steps per cycle, TT/closed form
is 0.974 at x = 100 and 0.977 at x = 150. Stronger absorption keeps ν smaller, so the
non-linearity is ~2 % there, against 25 % at 17.48 keV.

**Update 2 — a few-cycle packet in full TT (`explore/nb/e65_packet.py`, `e65_packet.log`).**
- Setup: the Entry 59 packet with α = 40, κ = 6 (g = 0.024), R_D = 0.2 (x ≈ 5), S = 0.02 / 0.2 / 0.7.
  The table gives W/S relative to the S = 0.02 run.

| S | ν estimate (πα\|J_1\|τ_d) | ν² | at η = q₀ (spectral peak) | at η = 2q₀ (twice the peak) |
|---|---|---|---|---|
| 0.2 | 0.19 | 0.035 | 0.973 | 1.221 |
| 0.7 | 0.62 | 0.385 | 0.738 | 2.604 |

- **Two different non-linearities.**
  - **(i) At the peak: saturation.** The weight falls by about ν² (2.7 % against 3.5 %; 26 % against
    38 %). This confirms the Entry 65 rule for packets, with ν = πα|J_1(S)|τ_d.
  - **(ii) In the wings: harmonic generation.** At 2q₀ the weight grows to 1.2× and 2.6× its linear
    value. This is the J₂(S) term of the Bessel expansion, a kinematic effect, already listed as
    candidate bias (ii) in Entry 57. It depends on S itself, not on ν. Where the linear spectrum is
    small, a harmonic of the strong peak dominates.
- **Chalise's case.** In physical units, ν ≈ π|J_1(S)|ℓ/Λ ≈ 0.03 (ℓ ≈ 0.3 µm), so (i) is
  negligible. But S ≈ 0.7 (Entry 57), so (ii) is not.
  - **Hypothesis, not checked:** the second harmonic of their 87 GHz peak sits at 174 GHz, which is
    where Entry 57 saw an unexplained feature in their reflectivity spectrum (150–170 GHz).
  - The difference frequency 113 − 87 = 26 GHz falls in their low-frequency band too.
  - Checking this needs their packet and energy, so it is not done here.
- **Rule for the estimator, both parts.**
  - Saturation is controlled by ν ≲ 0.2.
  - Harmonics are controlled by S itself (|J₂/J₁| ≈ S/4).
  - Only the first is a dynamical-diffraction effect. The second is kinematic and known in kind.

---

## 2026-10-15 — Entry 66: Chalise et al. Fig. 2 tails in full TT at Si 400 — second-order strain terms, not the crystal

**Why.** Entry 65 left a hypothesis unchecked: the 150–170 GHz feature in their Fig. 2b may be the second
harmonic of the 87 GHz packet. Entry 65 said this needed their packet and energy. A semi-infinite solver
removes the cost problem, so the case is run here at the real crystal (Si 400, 17.48 keV) with an analytic
packet of their peak strain. Script `explore/nb/e66_harmonic.py` (`check`, `run`, `fig`); data
`e66_harmonic.json`, log `e66_harmonic.log`; figure `figures/si400_harmonic.png`.

**1. Method: one ODE gives the whole R(t).**
- The packet enters at t = 0 and moves at v = 1, so phi(zeta, t) = S g(t − zeta) depends on depth only
  through s = t − zeta, the distance behind the front.
- The Bragg Riccati equation has no other explicit depth dependence. Below the front the crystal is perfect,
  so R = r0 there (Darwin–Prins root).
- Integrate dR/ds = −dR/dzeta from s = 0 upward. The value at s = t is the surface reflectivity at time t.
  This needs no slab thickness and no loop over time. One run (126 frequencies, 30k steps) takes about a
  second.
- **Check** against `solve_bragg` (α = 8, κ = 6, S = 0.3, R_D = 0.2): max |ΔR| = 2.9e-4 at |R| ≈ 0.24.
  This equals the slab's back-face term with no packet (2.9e-4, 2.2e-4, 1.2e-4 at the three η). So the two
  solvers agree up to the finite slab itself.
- **Convergence:** 64 and 128 steps per cycle give identical K to 4 digits.
- Units: length Λ = 11.54 µm (α = 1), time Λ/v_L = 1.368 ns, κ = 2πg = 0.063.
  - S is set by peak strain 1e-3: S = ε GΛ/q0, which gives 0.71 at 87 GHz and 0.55 at 113 GHz.
  - The packet is the Entry 58 shape: φ = S e^{−s/τ_d} sin(q0 s), so the strain starts with a step, as in
    their Fig. 2a. The decay is 0.9 cycles (their peak FWHM ≈ 30 GHz) or 2 cycles.
  - Window 140 ps (their Fig. 2a: entry at 60 ps, end at 200 ps) and 1400 ps.
- **Measure:** K(f) = (W(f; S)/S) / (W(f; S/100)/(S/100)), the non-linear over the linear weight at the
  same frequency.
  - The linear transfer function (Entry 59/65 H, the heterodyne factor, the window) cancels in K.
  - So K is the strain non-linearity alone.

**2. Results (140 ps window; 1400 ps changes K by ≤ 0.1).**

| packet | S | K(f0) | K(2f0) | K max (at) | K at 20 / 30 / 40 GHz |
|---|---|---|---|---|---|
| 87 GHz, 0.9 cyc | 0.71 | 1.009 | 1.17 | 1.28 (158 GHz) | 0.88 / 0.81 / 0.81 |
| 87 GHz, 2 cyc | 0.71 | 1.036 | 1.74 | 1.85 (170 GHz) | 0.65 / 0.62 / 0.69 |
| 113 GHz, 0.9 cyc | 0.55 | 1.003 | 1.11 | 1.22 (202 GHz) | 0.93 / 0.86 / 0.84 |
| 113 GHz, 2 cyc | 0.55 | 1.012 | 1.51 | 1.66 (216 GHz) | 0.75 / 0.68 / 0.71 |

- **Scaling** (87 GHz, 0.9 cycles, S = 0.71 / 0.36 / 0.18): (K − 1)/S is constant at 140–158 GHz
  (0.32–0.40) and at 200 GHz (−0.42 to −0.52).
  - The deviation is the second-order (J₂ and mixing) term interfering with the first-order one.
  - It is not ν² saturation: at the peak K = 1.01, as Entry 65 predicted (ν ≈ 0.03).
- **High-frequency wing.**
  - The weight rises above the linear value from about 1.3 f0 up to 2 f0, by 20–85 % at the maximum.
  - It falls below the linear value just past 2 f0 (K ≈ 0.5–0.7 at 200–230 GHz for the 87 GHz packet).
  - In their normalised units, the shoulder lies at 0.2–0.35 of the peak, where the strain is 0.1–0.2.
- **Low-frequency wing.**
  - The difference-frequency terms reduce the weight to K = 0.6–0.9 at 20–40 GHz.
  - The longer packet gives the larger reduction.

**3. Against their Fig. 2 (figure downloaded from the arXiv HTML; numbers read off the image).**
- **Location.**
  - For 87 GHz, the computed shoulder ends at 158–170 GHz. They see a step in the reflectivity curve at about
    165 GHz.
  - For 113 GHz, the computed shoulder is at 200–216 GHz. Their reflectivity curve is above the strain curve
    from about 155 GHz upward.
  - The location therefore agrees with the hypothesis of Entry 65.
- **Size and sign.**
  - Their excess at 170–240 GHz is about +0.02 of the peak.
  - Here the excess is +0.05 to +0.1 of the peak before 2 f0, and the deficit past 2 f0 is of the same order.
  - Their strain spectrum has a notch near 165 GHz (b) and 215 GHz (d), which the analytic packet does not
    have. The sign past 2 f0 depends on the packet's phase structure, so a shape match is not expected.
- **Low frequencies: two factors multiply.**
  - At 20 GHz, the crystal factor of Entry 65 gives H(20)/H(100) = 0.51.
  - The strain non-linearity gives K = 0.65–0.88.
  - Together that is 0.33–0.45, so the reflectivity under-reads 20 GHz by 2.2–3×. Their figure shows about 4×.
  - Both factors are in their forward model, since they solve TT with the full phase.

**Verdict.**
- At their strain (1e-3 at 87–113 GHz, S = 0.55–0.71), the reflectivity-spectrum estimator is linear at the
  peak (≤ 4 %), but not in the wings.
- Second-order strain terms change the weight by tens of percent at 1.3–2.5 f0 and at f ≲ 0.5 f0.
- The location of the extra weight agrees with the unexplained 150–170 GHz feature in their Fig. 2b.
- The hypothesis of Entry 65 moves from "unchecked" to "consistent in location and order of magnitude". The
  shape and sign are not tested, because their packet is not available.
- **Physics class.** K is kinematic in kind. It is the Bessel/mixing expansion, candidate bias (ii) of
  Entry 57, and is known in principle. What is new is the size at their stated strain and the separation of
  this effect from the crystal factor H.
- **Usable criterion (the Entry 57 stop rule).**
  - The non-linear part grows as S while the linear part does not.
  - Two runs or measurements at different fluence would separate the two parts.
  - Scaling (W/S at f) between fluences changes by (K − 1) ∝ S in the wings, and stays flat where the
    spectrum is linear.
  - This criterion needs only a fluence series, not their packet.

**Not done.**
- Their actual udkm1Dsim packet. The notch in their strain spectrum and the exact shoulder shape need it.
- How they combine angle-resolved FFTs into one spectrum (Entry 57). Here the weight is taken at ω = |η|,
  as in Entry 58.
- Instrument resolution (their Sec. IV).


---

## 2026-10-19 — Entry 67: amplitude–rise-time map of envelope-induced overshoot

**Question.** Where does a linearly increasing acoustic amplitude produce a transient
first-satellite centre intensity above the fully filled steady value?

**Setup.** `explore/nb/e67_overshoot_map.py`; Laue, alpha = L/Lambda_P = 0.5,
alpha_ac = L/Lambda_ac = 40 (R = 0.0125), n = 1, no absorption. The envelope
is a causal *linear ramp*, A(s) = clip(s/T, 0, 1), held at 1 afterwards.
This is the ramp used by `memory_drive.ramp`, not the smoothstep in `Drive`.
The scan covers S = 0.5–6 (111 points) and T/t_tr = 0.05–4 (101 points),
11,211 parameter pairs. The transient window is 0 <= t/t_tr <= 1 + T/t_tr;
the final uniform envelope fills the whole slab at its right endpoint.

**Observable and method.** Centre intensity at eta = q, normalized to incident
intensity, interpreted as the carrier-phase-averaged response. The map uses the
independent-satellite accumulated-coupling approximation from Entries 5 and 28:

    Theta_1(t) = integral_0^1 J_1(S A(t-zeta)) d zeta
    I_1(t) = sin^2(pi alpha Theta_1(t))
    I_steady = sin^2(pi alpha J_1(S))
    Delta I = max_t I_1(t) - I_steady

A cumulative Simpson antiderivative of J_1(S x), x in [0,1], evaluates the
sliding depth integral, including the causal zero region and the held final value.
Colour uses the absolute excess, not a ratio: the latter becomes ill-conditioned
near a dark steady satellite. This is a centre observable, not an integrated peak
population. At alpha = 0.5, |Theta_1| <= max|J_1| < 0.59, so the response does not
cross a Pendellösung maximum. The overshoot here comes from the envelope passing
through stronger Bessel coupling, not from a full intensity-exchange cycle.

**Results.** Data: `e67_overshoot.npz`; validation summary: `e67_overshoot.json`.
Poster figure: `figures/poster/f7_overshoot.png`, generated by `make_poster_figs.py f7`.

| case | S | T/t_tr | peak | steady | excess |
|---|---|---|---|---|---|
| A | 1.5 | 1 | 0.590505 | 0.590505 | 0 |
| B | 3.0 | 1 | 0.436665 | 0.257829 | 0.178836 |
| C | 3.85 | 2 | 0.532790 | 0.00013330 | 0.532657 |

B reaches 1.6936 times its final value. C is close to the first J_1 zero but
not exactly at it. The strongest excess on this sampled domain is 0.603419
at S = 3.8, T/t_tr = 4. This lies at the upper rise-time boundary and is not
an unconstrained optimum. Longer ramps allow more of the crystal to sample
amplitudes near the Bessel maximum simultaneously. Dependence on S is not monotonic.

**Checks.** The entire map was recomputed with 2,001 instead of 1,001 time
samples and 8,193 instead of 4,097 antiderivative samples. The largest change
in Delta I was 4.50e-6. The saved map is the refined one.

Full TT checks use S = {1.5, 3, 3.85, 5}, T/t_tr = {0.2, 1, 3}, 101 times
per trace, eight equally spaced carrier phases at each *fixed envelope time*,
and `required_nz` depth resolution. The largest absolute trace error against
the approximation was 8.65e-5; the largest error in excess was 8.18e-5.
For A–C, doubling both carrier phase samples (8 to 16) and depth resolution
changed phase-averaged TT intensity by at most 3.41e-5. The plotted TT dots
use the refined 16-phase calculation. These checks are selected conditions,
not a full TT sweep or a global error bound for the entire map.

**Audit.** `e67_worker_audit.py` (results `.json`, runs in `.log`) rechecks A–C,
the map maximum (S = 3.8, T/t_tr = 4) and the short-rise corner (S = 6,
T/t_tr = 0.05). Direct `quad` depth integration matches the antiderivative law
to 1.7e-8; causality and the filled endpoint hold exactly. Full TT with 16
phases: depth doubling changes intensity by at most 3.4e-5 (A), 16 to 32 phases
by at most 1.5e-8, and the TT minus approximation error is at most 8.3e-5.
Excess errors are at most 5.2e-5 (B). The 6.0e-4 excess at the short-rise corner is
comparable to this error scale and is not a map feature. The map maximum is
still rising at T/t_tr = 4 toward the quasi-static bound
sin^2(pi alpha max|J_1|) = 0.6272. B (T = 1) and C (T = 2) lie between map
rows (0.998, 1.9855). Their table values are computed at the exact T.

**Scope.** The map establishes behaviour for this fixed thickness, geometry,
order, and linear envelope. It does not map thickness-dependent Pendellösung
overshoot, arbitrary envelope shapes, absorption, or an experimental detection limit.
The mechanism is the one already documented in Entry 15; the new result is its
parameter map and explicit checks, not a claim of a newly discovered mechanism.

---

## 2026-10-22 — Entry 68: is the linearity limit of the Fourier-weight estimator set by one number? (falsification run for the paper's central claim)

**Why.** `notes/paper_candidates.md` proposes as the paper's central claim that a satellite is a
dynamical reflection with its own amplitude, and that this amplitude decides whether time-resolved
Bragg reflectivity is a linear spectrometer of the strain (Entry 59/65 transfer function H) or a
dynamically forming satellite (Entry 1/9 law). Isaac's question B3 asks whether the short-packet
linear regime and the long-train dynamical regime are joined by **one** dimensionless criterion.
Entry 65 proposed ν = πα|J_1(S)| ℓ/Λ on four points. Before building on it, try to break it.

**Hypothesis under test (strong form).** With ν defined as the first-order satellite amplitude at
the surface, ν = max_t |R_lin(t) − r0| (computed from a linear run and scaled by S/S_lin), the
saturation of the Fourier weight at the spectral peak, 1 − K with
K = (W(S)/S)/(W(S_lin)/S_lin), is a function of ν alone — the same for any packet length,
absorption and satellite position.

**Falsification rule, fixed before the run.** If at equal ν the value of 1 − K spreads by more than
a factor 2 across (g, x, ncyc), the strong form is false.

**Setup.** `explore/nb/e68_nu_collapse.py` (`run`, `stats`, `fit2`, `conv`, `window`, `fig`);
data `e68_nu_collapse.json`; logs `e68_nu_collapse.log`, `e68_stats.log`; figure
`figures/paper/e68_nu_collapse.png`. Reproduce: `uv run python explore/nb/e68_nu_collapse.py run stats fit2 fig`
(28 s for `run`).
- Semi-infinite symmetric Bragg, σ, normal absorption only, retarded-time ODE of Entry 66 (one
  integration gives R(t)); length unit Λ (α = 1), κ = 2πg.
- Packet φ = S g(s), g = e^{−s/τ_d} sin(q0 s), τ_d = ncyc/x, q0 = 2πx (satellite at x Darwin
  half-widths). Strain starts with a step at the front, as in Entries 58–66.
- Grid: g ∈ {0.01, 0.046, 0.12} (Si 400 at 17.5 keV, 8 keV, and Entry 59), x ∈ {2, 5},
  ncyc ∈ {1, 3, 10, 30}, S ∈ {0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.5, 1}, S_lin = 1e-3. 24 series,
  192 nonlinear points. Weight read at η = q0 (also 0.8 q0 and 1.25 q0, stored, not analysed).
- Window T = 5τ_d + 4/Re(Λ − iW) + 1. RK4 at 64 steps per cycle of (η + S q0).

**Results.**

1. Spread of 1 − K at fixed ν (log-log interpolation within each series, all 24 series):

| ν | min | max | max/min | least / most saturating series |
|---|---|---|---|---|
| 0.1 | 0.0012 | 0.0175 | 14.1 | g=0.046, x=2, 30 cyc / g=0.01, x=2, 1 cyc |
| 0.2 | 0.0050 | 0.0611 | 12.1 | same / g=0.01, x=5, 1 cyc |
| 0.3 | 0.0116 | 0.1116 | 9.6 | same / g=0.01, x=2, 1 cyc |
| 0.5 | 0.0370 | 0.2496 | 6.8 | same / g=0.01, x=5, 3 cyc |
| 0.8 | 0.1151 | 0.4227 | 3.7 | same / g=0.01, x=5, 10 cyc |

**The strong form is false.** At equal ν the saturation varies by an order of magnitude. Short
packets (which need a large S to reach a given ν) saturate most; the longest packets least.

2. It is not the local Bessel factor. Dividing K by 2J_1(S)/S widens the spread (26× at ν = 0.1),
so the extra saturation of short packets is not J_1 saturation of the drive.

3. Not a two-term law either. Fitting 1 − K = a S² + b ν² over the 144 points with S ≤ 0.3 gives
a = 0.113, b = 0.49 with pred/actual 5–95 % = [0.35, 4.0]; forcing b = 2/3 (the Bragg tanh² law)
drives a negative. A ν-only fit (b = 0.52) has 5–95 % = [0.19, 2.5]. The saturation depends on the
packet's shape relative to the extinction structure, not only on its peak amplitude and on S.

4. What survives, as a bound. Over the scanned domain, 1 − K ≤ 0.0175 at ν = 0.1, ≤ 0.061 at
ν = 0.2, ≤ 0.11 at ν = 0.3, i.e. **1 − K ≤ 1.8 ν² for ν ≤ 0.3** in every series. So "ν ≲ 0.1
guarantees linearity to 2 %" holds as a sufficient condition; ν is a usable safety number, not a
crossover law. In the figure the Bragg law 1 − tanh²ν/ν² runs through the middle of the cloud
within about a factor 3 for ν ≥ 0.1; the Laue sin² law is its lower edge.

5. Same ν, different route (g = 0.01, x = 5, S ≤ 0.3): at ν ≈ 0.3 the 1-cycle series is off the
grid (its S ≤ 0.3 points reach only ν = 0.09), the 10-cycle series gives 0.073, the 30-cycle series
0.053. The route matters.

**Checks.**
- Step convergence: 64 vs 128 steps per cycle at g = 0.01, x = 2, 30 cycles, S = 0.1 and 0.3: K
  differs by < 3e-7.
- Window: doubling T at g = 0.01, x = 5 changes K by ≤ 4e-4 (1 cycle) and ≤ 1e-4 (30 cycles).
- The retarded-time solver was checked against `solve_bragg` in Entry 66 (agreement to the slab's
  back-face term, 3e-4).
- The ν definition and the factor-2 rule were written in the script docstring before the run.

**Scope.** Bragg, semi-infinite, symmetric, σ, normal absorption only (no χ_ih), one packet family
(exponential decay, step strain front), S ≤ 1, satellite at 2 and 5 Darwin half-widths only (the
inside-Darwin regime x < 1 is not scanned), weight at the spectral peak only, ν taken as the
peak-in-time first-order amplitude. Not tried, deliberately (no second definition after the fact):
a time-averaged ν, other packet shapes, Laue geometry.

**Consequence for the paper (see `notes/paper_candidates.md`).** Question B3 gets a negative answer
in its strong form: there is no single dimensionless number that maps the linear transfer function
onto the dynamical satellite growth. The satellite's own amplitude ν is the right *quantity*
(it is exactly the formation law's phase πα|J_n|t̃ in Laue), it gives a sufficient linearity bound,
and it orders the regimes — but the saturation function is not universal. The central claim is
reworded accordingly: H is the correction, ν ≲ 0.1 is the guarantee, and the crossover is a region,
not a line.

---

## 2026-10-22 — Entry 69: Sondhauss & Wark (2003) vs Entry 59, aligned term by term (re-review, no recomputation)

**Why.** BACKLOG T1a asks for the Entry 64 comparison with normalisation, boundary conditions,
absorption, observable and approximations made explicit, and for one paragraph on what Entry 59
adds. Nothing is recomputed: the numbers are those of `explore/nb/e64_sw.log` (2026-10-14) and
Entry 65. **PDF not re-checked tonight**: `SondhaussWark2003_ActaA59_7.pdf` is untracked and absent
from this worktree and from the shared checkout, so the S&W column is Entry 64's transcription of
their eq. 25–31 and p. 12. Entry 64 §5 already found that the two-field S&W form with the Bloch
wave restored equals the Entry 59 Riccati form to 6e-16, so the only modelling difference between
the two first-order theories is the depth ansatz; the rest of the table is bookkeeping that a
reader must have to compare numbers.

**1. Alignment table.**

| item | Sondhauss & Wark 2003 (transcription, Entry 64) | Entry 59 / 64 (this work) | matched in e64 by |
|---|---|---|---|
| equation | time-dependent TT (their eq. 14/16), two-beam eq. 26 with explicit ∂/∂t terms | quasi-static two-beam TT, frozen φ per instant (DESIGN.md) | dropping their 2Ω/ω terms; justified by β = 2πfL/c ≤ 1.2e-2 (Entry 62). Their eq. 36 Doppler shift of the satellite is not represented |
| geometry, polarisation | InSb 111 symmetric Bragg, 10 µm, one polarisation | symmetric Bragg, σ, dimensionless α, κ; slab with back face hidden by e^−κ (κ = 6) or semi-infinite (Entry 66) | same dimensionless equation; InSb numbers not reproduced |
| drive | one plane-wave phonon, δχ_h cos(Ωt − q·r); δχ from the displacement to first order | φ = S cos(W(t − ζ)), single mode, phase only (no δχ_h change of magnitude) | identical at first order in u: δχ_h = iχ_h h·u |
| perturbation order | first: δχ·δD dropped (eq. 25) | first: Riccati linearised about R0 | same |
| depth ansatz of the perturbation | constant d_h^±, constant perfect-crystal (D_0^p, D_h^p) (eq. 28–29) | δ(ζ) = (A_s, B_s) e^{(λ − isW)ζ} with λ = cR0, i.e. the actual semi-infinite Bloch wave, and the surface boundary condition | the **only** physical difference; kinds SW / SWn / Bloch in `e64_sw_compare.py` |
| normalisation of the satellite amplitude | δD_h itself (eq. 31) | δR = B − R0 A, normalised to the incident wave | kind SWn removes this difference; it changes little (Entry 64 table) |
| boundary conditions | implicit (plane waves, no surface condition on δ) | D0(0) = 1, R(1) = 0 → Darwin–Prins root R0 for the unperturbed state; δ(1) = 0 | Bloch kind |
| absorption | complex χ_0 inside the deviation parameter β_h (implicit); χ_h'' present in their χ_h | normal absorption κ only (η → η + iκ); χ_ih **not** included | κ = 6, α = 8 (g = 0.12); not matched to InSb |
| deviation variable | β_h (eq. 27, 33) | η (Darwin half-width 2πα); x = η/2πα | linear map |
| observable | rocking curve of \|D_h^p + δD_h\|² at an instant (Fig. 3–4); no Fourier analysis | F(η) = mean_t[(\|R\|² − ⟨\|R\|²⟩) e^{iWt}] = conj(R0) δR_− + R0 conj(δR_+); H = \|F\| / kinematic law | F built from each kind's δR_± |
| stated validity | satellite ≤ 10 % of main peak; sidebands farther than ~1 rocking width from the main peak; crystal thick enough for no Laue fringes (p. 12) | ν = \|δR_lin\| ≲ 0.1 (Entry 68 bound; ≲ 0.2 for ≲ 4 %, Entry 65); any η including the plateau; slab fringes suppressed by κ | their 10 % rule ≈ ν² ≲ 0.1 × I_main, the same kind of condition |
| verification | numerics for InSb, Fig. 3 (agreement), Fig. 4 (failure near the main peak) | full TT, nt × 2 and n_z × 2 convergence 1.4e-5 of max \|F\| | `e64_sw.log` |
| intended use | least-squares fit of a sum of mode solutions to a rocking curve → mode occupations (one paragraph, no demonstration) | divide the measured Fourier weight by H(ω; g) → strain spectrum; linearity bound ν | — |

**2. The matched numbers (from `e64_sw.log`, α = 8, κ = 6, S = 0.02), ratio |F_pred|/|F_TT| at the satellite η = ±W:**

| W / Darwin half-width | 20 | 5 | 2 | 1 | 0.5 |
|---|---|---|---|---|---|
| SW (constant amplitudes, δD_h) | 1.022 | 1.220 | 1.687 | 3.014 | 8.45 |
| SWn (constant amplitudes, normalised) | 1.021 | 1.207 | 1.567 | 1.673 | 12.96 |
| Bloch = Entry 59 Riccati | 1.005 | 1.003 | 1.001 | 1.000 | 1.000 |
| SW worst near the main peak, \|η\| < 1.5 w_D | 3.5 | 3.7 | 4.4 | 6.2 | 9.4 |

The first-order theories agree with each other and with TT where S&W say theirs is safe (2 % at 20
half-widths). Inside about five half-widths they part, and the difference is the depth ansatz
alone: normalising (SWn) does not repair it.

**3. What Entry 59 adds, in one paragraph.** The linearisation of the Takagi–Taupin equations
about the perfect-crystal solution, the superposition of modes, and the idea of reading phonon
content from the linear solution are Sondhauss & Wark's (2003) and are cited as such. Entry 59
keeps the depth dependence of the perfect-crystal Bloch wave and the surface boundary condition in
the perturbation, which S&W's constant-amplitude ansatz drops; that single change makes the
first-order solution hold to 1 % at every detuning, including the Darwin plateau and the near-peak
region where S&W's own Fig. 4 and our e64 table show factors of 3–13. On top of that, Entry 59
changes the observable from a rocking curve at an instant to the time-Fourier component of the
reflectivity, and writes its magnitude as the strain spectrum times a crystal transfer function
H(ω/ω_D; κ/2πα) that depends on the crystal only: it collapses onto one curve in ω/ω_D for
different packets, depends strongly on absorption, and is a factor 2 across 20–100 GHz for Si 400
at 17.5 keV (Entry 65). Entry 68 adds the validity bound ν ≲ 0.1 for the estimator. None of this
touches the S&W time-derivative terms (Doppler shift, polaritons), which are outside the
quasi-static model.

**4. Still unresolved.** Whether an earlier Bloch-wave first-order form exists in Chukhovskii &
Förster 1995 or Wark & Lee 1999 (unread); anomalous absorption χ_ih; asymmetric geometry; the
Doppler term of their eq. 36 at GHz (expected ~Ω/ω ~ 1e-8, negligible, not computed).

---

## 2026-10-22 — Entry 70: signal maximum vs Fisher-information maximum, one case (directive B2)

**Why.** Isaac's question B2: is the detuning where the signal is largest the detuning where the
information about the acoustic parameters is largest? One case, as boss asked: the Entry 66 packet
(87 GHz, decay 0.9 cycles, peak strain 1e-3, S = 0.71) at Si 400, 17.48 keV (g = 0.0100), with 8 keV
(g = 0.0458) as the crystal-factor control. Script `explore/nb/e70_fisher.py` (`run`, `diag`, `fig`);
data `e70_fisher.json`, `e70_fisher_f30.json`; logs `e70_fisher.log`, `e70_fisher_f30.log`; figures
`figures/paper/e70_fisher.png`, `e70_fisher_f30.png`. Reproduce:
`uv run python explore/nb/e70_fisher.py run fig` (24 s) and `E70_FMIN=30 uv run python explore/nb/e70_fisher.py run fig`.

**Measurement model.** R(t) = |R|² at one detuning η (one angle), sampled uniformly over a window
T (140 ps as in their Fig. 2a; 1400 ps check), additive Gaussian noise of constant σ. Signal:
W(η) = Fourier weight at ω = |η| (Entries 58–66). Information: per-detuning Fisher matrix
I_ij(η) = σ⁻² Σ_t ∂R/∂θ_i ∂R/∂θ_j for θ = (S, f0, τ_d), central differences (1 %; 0.5 % changes
I_ii by ≤ 4.6e-3). Reported: argmax of I_ii and of the marginal precision 1/[I⁻¹]_ii.

**Pre-registered rule (in the script docstring before the run).** For S, "different" only if
|f_I(S) − f_W| exceeds both one Darwin half-width (0.73 GHz) and 10 % of f0 (8.7 GHz). For f0 and
τ_d a flank maximum is expected from Fourier theory and is not a finding; a crystal effect is
claimed only if the 17.48 keV maximum differs from the 8 keV one by more than 10 % of f0.

**Result on the pre-registered grid (10–260 GHz).** W and I_SS both peak at the grid edge, 10 GHz,
at both energies and both windows. `diag`: the linear run (S/100, scaled) gives the same edge
maximum (W(10)/W(87) = 3.05 vs 2.92), so this is not the second-order S² term; it is the
first-order low-frequency rise of W ∝ |r0||φ̂| ∝ |ε̂|/ω² (Entry 59 §1, strain step at the front)
read through a window holding 1.4 cycles at 10 GHz. **Post-hoc change, labelled as such:** the grid
is cut at 30 GHz (≥ 4.2 cycles in 140 ps). Everything below is from that rerun.

| | 17.48 keV, 140 ps | 17.48 keV, 1400 ps | 8 keV, 140 ps | 8 keV, 1400 ps |
|---|---|---|---|---|
| f_W (signal max) | 84 | 84 | 84 | 84 |
| f_I(S), argmax I_SS | 87 | 87 | 87 | 87 |
| f_I(f0) | 79 | 80 | 79 | 80 |
| f_I(τ_d) | 87 | 87 | 87 | 87 |
| marginal precision S / f0 / τ_d | 30 / 87 / 85 | 30 / 88 / 49 | 30 / 87 / 85 | 30 / 88 / 49 |

**Reading.**
- Amplitude: f_I(S) − f_W = 3 GHz, below the 8.7 GHz threshold. **They coincide.** B2 is answered
  "no" for the amplitude in this case, as the trivial expectation said (R is linear in S).
- Frequency: the information sits on the low flank (79–80 GHz against a peak at 84), which is
  Fourier theory (the f0 derivative of a peaked spectrum is largest on the flanks); the marginal
  precision, which has to separate f0 from S and τ_d, moves back to the peak (87–88). Not a finding.
- Crystal: the 8 keV and 17.48 keV columns are identical to the GHz. **No crystal effect** on any
  maximum, far below the 10 % criterion. Entry 71 explains why: in 140 ps to 1.4 ns the packet is
  0.1 to 1 Λ deep and neither absorption nor the extinction phase has acted.
- The marginal-precision maxima at the grid edge (S at 30 GHz, τ_d at 49 GHz for 1400 ps) are
  where the three derivatives are least collinear, not a signal feature; not interpreted.
- The kinematic H = 1 reference (|r0| |φ̂| and its derivatives) also peaks at the grid edge after the
  φ̂ correction, so it does not discriminate here; it is kept in the JSON, not used.

**Verdict.** For a picosecond packet read at one angle, the signal maximum and the Fisher maximum
for the amplitude coincide, the frequency information is on the flank as in any spectrum, and the
crystal moves nothing. **B2 yields no paper figure from this case.** The place where B2 could still
bite is a long train at a weakly absorbing reflection (window ≫ 1/κ), where H(ω) is not flat
(Entry 71); not computed.

**Scope.** Gaussian noise, constant σ, one angle at a time (no angle-resolved joint FIM), one packet
shape, S = 0.71 (second-order wings present, Entry 66), semi-infinite symmetric Bragg, σ, normal
absorption only.

---

## 2026-10-22 — Entry 71: the crystal factor H needs a long window; Entry 65 §2 corrected for the Chalise case

**Why.** Entry 70 found the 8 keV and 17.48 keV responses of a 140 ps packet measurement identical,
although Entry 65 gives H(20 GHz)/H(100 GHz) = 0.51 at 17.48 keV and 0.93 at 8 keV. Entry 65 computed
H for a steady single mode filling the crystal; Entry 65 §2 and Entry 66 §3 then applied that number
to Chalise et al.'s 140 ps measurement. Test whether that step was allowed. Script
`explore/nb/e71_window_H.py`, data `e71_window_H.json`, log `e71_window_H.log`. Reproduce:
`uv run python explore/nb/e71_window_H.py` (11 s).

**Setup.** Full semi-infinite TT (Entry 66 solver), the Entry 66 packet at S/100 (linear), windows
140 ps, 1.4 ns, 14 ns; H_meas(f) = W(f) / (|r0(f)| |φ̂(f)|) with φ̂ = ε̂/ω (Entry 59 §1), normalised
at 100 GHz; compared with the steady closed form H(f)/H(100) of Entry 65.
**Pre-registered rule:** Entry 65's number for their case stands only if H_meas(20)/H_meas(100) at
140 ps is within 20 % of 0.51 (17.48 keV); if it is near 1, Entry 65 §2 is wrong for short windows.
(The first run of the script normalised by |ε̂| instead of |φ̂|; the printed ratios were the correct
ones times 100 GHz/f. Fixed before the numbers below; the rule was not changed.)

**Result**, H_meas(f)/H_meas(100 GHz):

| f GHz | 20 | 30 | 40 | 60 | 87 | 150 |
|---|---|---|---|---|---|---|
| steady H, 17.48 keV (Entry 65) | 0.511 | 0.676 | 0.786 | 0.909 | 0.982 | 1.034 |
| 140 ps window, 17.48 keV | 1.015 | 1.043 | 1.039 | 1.009 | 0.957 | 1.039 |
| 1.4 ns | 0.999 | 0.999 | 1.000 | 0.999 | 0.996 | 1.006 |
| 14 ns | 0.944 | 0.975 | 0.986 | 0.994 | 0.998 | 1.004 |
| steady H, 8 keV | 0.929 | 0.968 | 0.983 | 0.994 | 0.999 | 1.002 |
| 8 keV at 20 GHz, 140 ps / 1.4 ns / 14 ns | 1.014 / 0.998 / 0.958 | | | | | |

**The rule fires: Entry 65 §2 is wrong for a 140 ps measurement.** The crystal factor is flat to
±4 % (140 ps, leakage) and to 0.1 % (1.4 ns) at both energies. The steady suppression begins to
appear only at 14 ns (0.944 against 0.511 steady).

**Why, in one line.** The steady H = κ/|κ − 2cR0| with |2cR0| ≈ π/x (per Λ, x = f/f_D) is the
modulus of the infinite depth integral of e^{−(κ + iπ/x)ζ}. Over a finite depth T (the window in
units of Λ/v_L = 1.368 ns) the integral is (1 − e^{−(κ+iπ/x)T})/(κ + iπ/x), which tends to T for
small T: no suppression. Hand check at T = 10.3 (14 ns): the ratio of this expression at 20 GHz to
that at 100 GHz is 0.947, against 0.944 from TT. So H needs the window, and hence the strained
depth, to exceed both 1/κ (16 Λ = 185 µm = 22 ns at 17.48 keV) and x/π (8.7 Λ at 20 GHz). A
picosecond packet 0.1 Λ deep sees none of it.

**Consequences.**
- Entry 65 §2 ("a peak-normalised Fourier-weight spectrum taken at 17 keV under-reads 20 GHz by a
  factor 2") holds for a **steady or long acoustic train** read over tens of ns, the poster's own
  GHz-transducer geometry, and **not** for the Chalise et al. 140 ps packet. Entry 66 §3's "two
  factors multiply, 0.33–0.45" is withdrawn: only the second-order factor K = 0.65–0.88 applies
  to their window, which leaves most of their ~4× low-frequency deficit unexplained by this model
  (instrument, their windowing and normalisation are the remaining candidates, Entry 57).
- `notes/qa_prep.md` §6b and the poster Outlook sentence quoting the Si 400 numbers are right only
  with the window condition attached. The poster is frozen; the correction goes to
  `notes/paper_candidates.md` and the paper.
- Entry 59's closed form itself is not wrong: it is the infinite-window limit, and its TT checks
  (κ = 6, T ≥ 0.6 L ≫ 1/κ) satisfied the condition. The window condition was never stated; it is now.
- For the paper this moves H from "correction for picosecond-ultrasonics spectra" to "property of
  long-train (GHz) measurements", i.e. into the same experiment as the formation law. The two
  halves of P1/P2 now share one window condition and one apparatus.

**Not done.** A windowed closed form (the finite-T integral above, with G(W) for the packet);
the intermediate regime T ~ 1/κ; windows beyond 14 ns (the packet is then 10 Λ deep and the
quasi-static and plane-wave assumptions should be rechecked).

---

## 2026-10-22 — Entry 72: which histories time-resolved diffraction does not distinguish (directive B1)

**Why.** Isaac's question B1: can two acoustic fields that steady-state diffraction cannot tell apart
be told apart in time-resolved diffraction? Inside this model the steady state of a filled
constant-amplitude wave knows only (S, q); any two drive histories with the same final amplitude are
indistinguishable there, and Entries 28–31 already show them distinguishable in the transient. The
quantitative content of B1 is therefore the converse: **what does the transient not distinguish?**
By the accumulated-coupling law, order n reads only the running mean of J_n(S A(s)) over the last
transit. Script `explore/nb/e72_history_degeneracy.py` (`law`, `tt`, `fig`); data
`e72_history_degeneracy.json`; log `e72_history_degeneracy.log`; figure
`figures/paper/e72_history_degeneracy.png`. Reproduce: `uv run python explore/nb/e72_history_degeneracy.py` (7 s).

**Construction.** History 1: linear ramp A_1(s) = clip(s, 0, 1) (T = t_tr, then held). History 2:
J_1(S A_2(s)) = J_1(S A_1(s)) + δ(s), δ = 0.1 sin(2π(s − 0.2)/0.5) on s ∈ [0.2, 0.7], zero outside;
inverted on the monotone branch of J_1 (S = 1.5, S A ≤ 1.84). The envelopes differ by up to 0.167
(a 17 % amplitude wobble during the ramp; A_2 is non-monotone). For t ∈ [0.7, 1.2] the window
[t − 1, t] holds the whole wobble, whose mean is zero, so Θ_1 is identical there by construction.
Both histories end in the same filled state S = 1.5.

**Pre-registered rule (docstring, before the TT run).** Laue, α = 0.5, α_ac = 40 (R = 0.0125,
Entry 67 conditions), 8 carrier phases at fixed envelope time, η = +q and +2q, n_z = 2800:
(i) |ΔI_1| < 1e-3 for t ∈ [0.7, 1.2] → order-1 degeneracy confirmed, else claim #21 fails here;
(ii) TT difference traces agree with the law's to 1e-3 at every t, both orders;
(iii) if the law predicts max|ΔI_2| > 3e-3 in the window and TT confirms, order 2 breaks the
degeneracy at this S; if the law predicts less, "not resolved at this S", no claim.

**Results.**

| | law | full TT |
|---|---|---|
| max\|ΔI_1\|, t ∈ [0.7, 1.2] | 8.4e-11 | **8.8e-7** |
| max\|ΔI_1\|, outside the window | 2.50e-2 | 2.50e-2 |
| max\|ΔI_2\|, in window | 1.25e-3 | 1.25e-3 |
| max\|ΔI_2\|, outside | 5.82e-3 | 5.82e-3 |
| max\|ΔI_TT − ΔI_law\|, n = 1 / 2 | — | 2.5e-6 / 4.4e-6 |
| max\|I_TT − I_law\| (each history), n = 1 / 2 | — | 1.9e-5 / 6.7e-5 |
| convergence (8 → 16 phases, n_z × 2), n = 1 / 2 | — | 5.7e-6 / 5.3e-5 |

- (i) **Confirmed.** Two histories whose envelopes differ by 17 % give the same first-order
  satellite to 9e-7 while the window contains the perturbation, and differ by up to 2.5e-2 before
  and after. The order-1 transient reads only the t_tr running mean of J_1(S A), as the law says.
- (ii) **Confirmed** at 4e-6, three orders below the threshold, on a non-monotone envelope.
  Claim #21 (Entry 33) survives this attempt.
- (iii) **Not resolved at this S.** The law predicts 1.25e-3 for the order-2 difference inside the
  window, below the pre-set 3e-3; TT reproduces the number but no claim is made that order 2
  breaks the degeneracy. (At S = 1.5 the n = 2 satellite is weak, I_2 ≤ 0.085; a larger S is outside
  the monotone-inversion construction used here and was deliberately not tried tonight.)

**Answer to B1, bounded.** Yes, the transient distinguishes histories the steady state cannot, but
only up to a known equivalence class: at each order n it measures the moving average of J_n(S A)
over one transit, so histories that differ by a zero-mean (per window) perturbation of J_n(S A)
are degenerate at that order. Different orders have different kernels, so multi-order readout
reduces, and may remove, the degeneracy; the size of that effect at S = 1.5 is 1e-3 and was not
claimed. This is a consequence of the Θ_n law (confirmed phenomenon within the model), not a new
mechanism, and it is the precise sense in which "time-resolved diffraction reads the drive history":
through a t_tr low-pass filter, order by order.

**Scope.** Laue, no absorption, α = 0.5 (no Pendellösung crossing), R = 0.0125, one perturbation
shape, n = 1, 2 only. The Bragg reading tanh²(πα|Θ_n|) is exact at resonance (Entry 30), so the
same degeneracy holds there; not run.

---

## 2026-10-23 — Entry 73: what S = 0.5–6 means physically (T1c, directive C)

**Why.** The poster's maps run S = 0.5–6 without saying what displacement, strain, stress and drive
power that is. Script `explore/nb/e73_realisability.py` (constants only, no TT; 1 s), data
`e73_realisability.json`, log `e73_realisability.log`. X-ray constants from X0h (Entry 20);
c44 = 79.6 GPa, ρ = 2329 kg m⁻³, v_T = 5846 m/s (textbook, not re-verified); acoustic nonlinearity
parameter β ≈ 5 (order of magnitude, not verified). Geometry: symmetric Laue, bulk shear wave along the
surface normal polarised along h (Entry 6), L = 3 Λ_P.

**Definitions.** u0 = S d/(2π); strain ε = K u0 = S d/Λ_ac (Λ_ac = v_T/f); stress σ = c44 ε; energy
density w = ½ c44 ε²; intensity I = v_T w; second-harmonic strain after L, ε2/ε = (β/2) k ε L;
R = Λ_ac/Λ_P (independent-satellite rule R ≤ 0.1, Entries 19, 54); β_qs = 2πfL/c (Entry 62).

**1. Displacement vs strain.** [Corrected in Entry 80: u0 = 15.3 pm × S; the 7.6 is the S = 0.5 row.] u0 does not depend on frequency: u0 = 7.6 pm × S for Si 440 (d = 96 pm),
so S = 6 is u0 = 92 pm = 0.95 d — the atoms move by almost one lattice spacing of the reflection. The
strain is what elasticity sees, and it is small: ε = 1.6e-5 × S at 1 GHz, 3.3e-5 × S at 2 GHz (Si 440).
A large S is a large absolute displacement at a small strain, because Λ_ac/d ≈ 6×10⁴.

**2. Table, Si 440 (same for 17.48 and 30 keV except R, β_qs and ε2/ε, which scale with L).** [The ε2/ε columns and the harmonic parts of the verdict are withdrawn for this shear geometry in Entry 80.]

| f GHz | S | u0 pm | u0/d | ε | σ MPa | w J/m³ | I W/cm² | ε2/ε (17.48 keV, L = 163 µm) | ε2/ε (30 keV, L = 296 µm) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.5 | 7.6 | 0.08 | 8.2e-6 | 0.65 | 2.7 | 1.6 | 0.004 | 0.007 |
| 1 | 1 | 15.3 | 0.16 | 1.6e-5 | 1.3 | 11 | 6.3 | 0.007 | 0.013 |
| 1 | 2 | 30.6 | 0.32 | 3.3e-5 | 2.6 | 43 | 25 | 0.014 | 0.026 |
| 1 | 3 | 45.8 | 0.48 | 4.9e-5 | 3.9 | 97 | 56 | 0.022 | 0.039 |
| 1 | 4 | 61.1 | 0.64 | 6.6e-5 | 5.2 | 172 | 100 | 0.029 | 0.052 |
| 1 | 6 | 91.7 | 0.95 | 9.9e-5 | 7.8 | 386 | 226 | 0.043 | 0.078 |
| 2 | 1 | 15.3 | 0.16 | 3.3e-5 | 2.6 | 43 | 25 | 0.029 | 0.052 |
| 2 | 2 | 30.6 | 0.32 | 6.6e-5 | 5.2 | 172 | 100 | 0.058 | 0.105 |
| 2 | 3 | 45.8 | 0.48 | 9.9e-5 | 7.8 | 386 | 226 | 0.086 | 0.157 |
| 2 | 6 | 91.7 | 0.95 | 2.0e-4 | 15.7 | 1546 | 904 | 0.173 | 0.314 |

Other reflections at the same S and f: Si 400 has 1.41× and Si 220 2× the strain, stress and ε2/ε
(d larger), and 4× the intensity for 220. R at 1 / 2 GHz: 0.108 / 0.054 (440, 17.48 keV),
0.059 / 0.030 (440, 30 keV), 0.137 / 0.069 (400), 0.161 / 0.081 (220). β_qs ≤ 1.2e-2 everywhere.

**3. Sharp front vs smooth waveform.** A resonant transducer of quality factor Q rings up over ~Q/f.
At Si 440: rise/t_tr = 0.36 / 1.08 / 3.59 for Q = 10 / 30 / 100 at 1 GHz and 17.48 keV (0.18 / 0.54 /
1.79 at 2 GHz; 0.20 / 0.59 / 1.97 at 1 GHz and 30 keV). The step-onset law (poster Fig. 1 collapse)
needs rise ≲ 0.2 t_tr, i.e. Q ≲ 5 at 1 GHz: unlikely for a resonant device. **The regime a real
transducer gives is the ramp regime of Entry 67 (T/t_tr = 0.2–4)**, already mapped by the Θ_n law
with ≤ 8e-5 error against TT, so no new TT run is needed to replace the sharp front by a smooth one.
Within a cycle, the waveform stays sinusoidal while ε2/ε is small; above ~10 % (S ≥ 4 at 2 GHz, or
S = 6 at 30 keV) the Entry 35 harmonic non-ideality applies (law error 0.014 at ε2 = 0.2).

**4. Two validities, separately.**
- **Elastic wave.** Linear elasticity holds throughout (ε ≤ 2e-4, σ ≤ 31 MPa, orders of magnitude
  below fracture); the first correction is finite-amplitude harmonic generation, ε2/ε ≤ 5 % for
  S ≤ 3 at 1 GHz and ≤ 2 at 2 GHz (17.48 keV), reaching 17–31 % at S = 6, 2 GHz. Attenuation at
  1–2 GHz over 0.3 mm is ≤ 2 % amplitude (Entry 35 estimate, textbook order). So the elastic side is
  **valid for all S ≤ 6 at 1 GHz, and for S ≤ 2–3 at 2 GHz as a clean sinusoid**.
- **TT equation.** Two-beam, symmetric, quasi-static (β_qs ≤ 1.2e-2), and the displacement enters
  the phase exactly for any u0 (the model never linearises in u). The independent-satellite reading
  needs R ≤ 0.1: 1 GHz at 17.48 keV sits at the edge (0.108); 2 GHz, or 30 keV, is inside. Not
  assessed: X-ray beam divergence against the Darwin width, lateral acoustic inhomogeneity.
- **Drive power.** The one number this model cannot judge. S = 1 at 1 GHz is 6 W/cm² of acoustic
  intensity, S = 6 at 2 GHz is 900 W/cm². A pulsed drive of ~100 ns (≥ 2 t_tr) at 1 kHz is a duty
  cycle of 1e-4, so average powers are mW/cm² and heating is not the issue; the peak conversion of a
  GHz bulk transducer is. **No verified transducer number is in hand**; this is a library/Isaac item.

**Verdict (Si 440, 17.48–30 keV).**

| S | 1 GHz | 2 GHz |
|---|---|---|
| 0.5–2 | **feasible**: ε ≤ 3.3e-5, I ≤ 25 W/cm² (pulsed), ε2/ε ≤ 3 % | **feasible**: ε ≤ 6.6e-5, I ≤ 100 W/cm², ε2/ε ≤ 10 % |
| 3–4 | **limit**: I = 56–100 W/cm², ε2/ε 2–5 % | **limit**: I = 230–400 W/cm², ε2/ε 9–21 % (harmonics must be modelled) |
| 6 | **limit**: I = 230 W/cm², ε2/ε 4–8 % | **not as a clean wave**: I = 900 W/cm², ε2/ε 17–31 % |

"feasible / limit" here is on the elastic and TT side with the drive power unverified; a transducer
delivering ≥ 50 W/cm² peak at 1–2 GHz would make the limit column feasible.

**Consequence for the paper.** The poster's S = 3–6 overshoot cases (Entry 67 B, C) are at the
limit at 1 GHz and need a transducer number; the S ≤ 2 regime (where the first satellite forms with
τ_1 = 17–30 ns) is where the proposal should be written. The ramp regime, not the step, is the
observable one.

---

## 2026-10-23 — Entry 74: directive questions B3 and B4, classified against what exists (no new computation)

**Why.** BACKLOG T1b asks for one entry per directive question with the result classed as
confirmed / in prior art / unverified. B1 is Entry 72, B2 is Entry 70. B3 and B4 are answered by
entries already in this notebook; per the handoff rule (do not re-open what is already tested, do
not recompute), this entry only states the classification and the evidence. Nothing was run.

**B3. Can the short-pulse linear transfer function and the long-train dynamical satellite growth be
joined by one dimensionless criterion?**
- Answer: **no in the strong form, yes as a bound.** Entry 68 (pre-registered, 24 series, 192 points):
  at fixed first-order satellite amplitude ν the saturation of the Fourier weight spreads by ×14
  (ν = 0.1) to ×3.7 (ν = 0.8) across packet length, absorption and satellite position; a two-term
  law a S² + b ν² does not capture it either. What holds is 1 − K ≤ 1.8 ν² for ν ≤ 0.3 on the scanned
  domain, so ν ≲ 0.1 is a sufficient linearity condition.
- Entry 71 adds the other half of the junction: the crystal factor H of the linear regime needs the
  window to exceed max(1/κ, x/π) Λ/v_s, which is tens of ns, i.e. the long-train regime. So the two
  regimes meet in one experiment (GHz train, tens of ns), not in one number.
- Class: the ν bound is **confirmed** (within the model, Bragg, one packet family); the universal
  crossover is **falsified**; the window condition is **confirmed** (Entry 71, hand check 0.3 %).
- Prior art: none found for the bound or the window condition; S&W 2003 give a 10 % satellite rule
  of the same kind (Entry 64 §1).

**B4. When the independent-satellite approximation fails, are there unexpected quantitative features?**
- Answer: **yes, one, already found and bounded: Bessel-zero resurrection.** Entries 45, 51, 52, 54:
  a satellite dark at a Bessel zero reappears through a third-order path with
  g_eff = πα (R/2)² C_n(S), population P = sin²[πα (R/2)² C_n(S)] ∝ α⁶/α_ac⁴ (confirmed without a fit
  to < 1 % for R ≤ 0.1, four decades; renormalised zero S* − j = (R/2)² C_n/J'_{−n} to < 0.3 % at
  R ≤ 0.67), derived independently as an effective 2×2 Hamiltonian (Entry 52) and verified against the
  exact Floquet solver (Entry 42).
- Also bounded: the approximation holds (population error < 5 %) for R ≤ 0.34–0.57 depending on S
  (Entry 54); the first corrections are inter-order coherence and forward-beam sidebands at O(R²)
  (Entry 51); no order-space localisation (Entry 51 C); the R ~ 1 main-line dark window is the known
  X-ray/neutron–acoustic resonance (Entin; Iolin–Entin; Raitman; Entry 50).
- Observability (Entry 54): the resurrected order dominates the phase-averaged exit intensity only
  within |α_ac − integer| ≲ 1.3e-3 and an angular half-width w/q ≲ 2e-4; its share of the
  |η − q| < 0.45 q integral is 4.4e-5. **Theory, not an experiment proposal** (Entry 51 E).
- Class: **confirmed within the model** (full TT + Floquet, no fit); **"possibly new" only** until
  Entin (JETP Lett. 26, 392 (1977); JETP 77, 214 (1979)) and Iolin–Entin (JETP 85, 1692 (1983)) are
  read in full (Entry 50, poster_outline "Before printing"). Not re-opened tonight.
- For the paper: P4 in `notes/paper_candidates.md` (N 4, E 1, G 2): a future project, not the thesis
  paper.

**What T1b leaves open.** B1's order-2 degeneracy breaking at larger S (Entry 72 iii); B2 for a
long train at weak absorption (Entry 70 verdict); B3's windowed closed form (Entry 71 "Not done");
B4's two unread Soviet papers.

---

## 2026-10-23 — Entry 75: the crystal factor as a closed form in the window length (T2a)

**Why.** Entry 71 showed that H of Entry 65 is the infinite-window limit and gave a hand estimate
for a 14 ns window (0.947 against TT 0.944). BACKLOG T2a: write that estimate as a closed form
H(f; T_w), find the window where H is flat and the window where the Entry 65 factor 2 is back, and
test the closed form against TT over the whole range. Script `explore/nb/e75_window_cf.py`.

**Pre-registered rule (written before the run).** The closed form is the Entry 59/65 first-order
form with each depth integral truncated at the window, i.e. 1/(Λ − isW) replaced by
[1 − e^{(Λ − isW)T}]/(Λ − isW) (derivation in §1 below; it assumes a packet short against
min(1/κ, x/π) and fully inside the window, the "CF-A" form; "CF-B" keeps the packet's Laplace
transform and has no such assumption). Comparison: full semi-infinite TT (Entry 66 solver) at
S/100 for the Entry 66 packet (87 GHz, 0.9 cycles), Si 400 at 17.48 keV and 8 keV, satellite
frequencies 20–150 GHz, windows 0.14 ns to 300 ns. The weight is the Fourier component at |η| of
|R(t)|² − |R0|² over [0, T_w] (static level subtracted, not the window mean, so that no leakage
term enters either side; the difference from Entry 71's definition is reported, not used in the
rule). Convergence checks: step 64 → 128 per cycle at one window, S/100 → S/1000 at one window;
the analytic continuation of the unstrained Riccati past the packet is checked against the RK4
integration at 14 ns (must reproduce Entry 71's 0.944 to 0.1 %).
- If |W_TT / W_CF-A − 1| > 1 % at any (f, T_w) with the convergence checks passing, CF-A is
  **rejected** as the closed form in T_w; CF-B is then tested with the same rule. If both fail,
  "no closed form in T_w at first order" is recorded and the rule is not changed.
- The flat window and the revival window are read from the closed form only after it passes:
  T_flat = largest window with H(20)/H(100) ≥ 0.9, T_half where the ratio is halfway
  (0.756 at 17.48 keV), T_90 where it is within 10 % of the steady 0.511. Each is confirmed by one
  TT point.

**Setup (as run).** Si 400 symmetric Bragg, 17.48 keV (κ = 0.0629 per Λ, Λ = 11.54 µm, Λ/v_L =
1.369 ns) and 8 keV; Entry 66 packet g = e^{−s/τ} sin(q0 s), f0 = 87 GHz, 0.9 cycles (τ = 10 ps),
plus 10 and 100 cycles (τ = 115 ps, 1.15 ns) for the scope of the packet-free form; peak strain
1e-5 (S/100 of Entry 66) and 1e-6 for the convergence check; f = 20–150 GHz; 13 windows 0.14–300 ns.
TT: the Entry 66 semi-infinite Riccati in s = t − ζ, RK4 (64 steps per cycle of the highest
frequency) through the packet to s_p = 30 τ; beyond s_p the crystal is unstrained and the Riccati has
constant coefficients, so u = (R − r1)/(R − r2) = u_p e^{Λ(s − s_p)} exactly, |R|² − |R0|² is a sum
of exponentials and its Fourier integral over [s_p, T_w] is analytic: every window costs nothing
beyond one RK4 through the packet (0.2 s; 17 s for the 100-cycle packet). Reproduce:
`uv run python explore/nb/e75_window_cf.py run posthoc fig` (110 s), data `e75_window_cf.json`,
log `e75_window_cf.log`, figure `figures/paper/e75_window_H.png`.

**1. The closed form.** With δ(t) = −ic(1 − R0²) S ∫₀ᵗ e^{Λ(t−s)} g(s) ds (Entry 59 linearisation,
Λ = i(η + iκ) + 2cR0), one integration by parts gives, for the window [0, T_w],

    F_T(W) = ∫₀ᵀ δ e^{iWt} dt = ic(1 − R0²) S [G_T(W) − e^{(Λ+iW)T} L_T(Λ)] / (Λ + iW),
    G_T(W) = ∫₀ᵀ g e^{iWs} ds,   L_T(Λ) = ∫₀ᵀ g e^{−Λs} ds,
    Ĩ(η)  = conj(R0) F_T(η) + R0 conj(F_T(−η))     (Fourier component of |R|² − |R0|² at η).

CF-B keeps G_T and L_T (analytic for the damped sinusoid). CF-A, the packet-free crystal factor, sets
L_T = G(−η), G_T = G (packet short against the window and against 1/|Λ − iη|), which is the Entry 65
cos-mode form with every kernel 1/(Λ − isW) replaced by [1 − e^{(Λ − isW)T}]/(Λ − isW); its
resonant term alone is Entry 71's hand estimate, H(f; T_w) ≈ H_∞(f) |1 − e^{−(κ + iπ/x) T_w}|.
H is the weight divided by |G_T(η)| |R0| πS/κ (kinematic law), so TT, CF-A and CF-B are compared as
absolute numbers, not normalised at 100 GHz.

**2. Checks.** Analytic continuation vs plain RK4 to 14 ns: ratio 1.00000–1.00049. Mean-subtracted
weight at 14 ns normalised at 100 GHz: 20 GHz 0.944, i.e. Entry 71's TT number reproduced. Step
64 → 128 per cycle: ≤ 2e-5. S/100 → S/1000: H changes by −0.19 % (20 GHz) to +0.27 % (150 GHz),
so TT at S/100 carries a 0.2–0.3 % second-order residual; the 0.78–0.80 % of CF-B for the 10- and
100-cycle packets falls to 0.08 % and 0.21 % at S/1000 (post-hoc run), so it is the same residual.

**3. Result: the rule fires on CF-A, CF-B passes.** Worst |H_TT/H_CF − 1| over all (f, T_w):

| packet | energy | CF-A (packet-free) | CF-B (packet kept) |
|---|---|---|---|
| 0.9 cycles, τ = 10 ps | 17.48 keV | **8.3 %** (0.14 ns, 87 GHz); 2.1 % at 0.5 ns; ≤ 0.6 % for T_w ≥ 1.4 ns | 0.30 % |
| 0.9 cycles | 8 keV | 8.3 %; same pattern | 0.30 % |
| 10 cycles, τ = 115 ps | both | 41 % | 0.78 % (0.08 % at S/1000) |
| 100 cycles, τ = 1.15 ns | both | 49 % | 0.80 % (0.21 % at S/1000) |

So **CF-A is rejected** as a 1 % closed form in T_w, per the rule; **CF-B is the closed form in
T_w**: exact to first order in S at every window from 0.14 ns to 300 ns, both energies, three packet
lengths, with the residual equal to the second-order term of TT. CF-A's failure is at windows
T_w ≲ 1 ns ≈ 100 τ only, and its cause is in the dropped term: L_T(Λ) − G(−η) ≈ i(Λ − iη) ∂G/∂W, i.e.
the packet's group delay t_g = ∂ arg G/∂W shortens the effective window to T_w − t_g. Post-hoc check
(written after the run, not pre-registered): CF-A with T_w → T_w − t_g(η) is within 1.3 % at 0.14 ns
and 0.5 % from 0.5 ns on, so the mechanism is identified but the Taylor form is not a 1 % closed form
at 0.14 ns either (the remaining part is the phase of the −η image term, which CF-A takes from a
cos mode). Nothing was changed in the rule; the packet-free form stands only for T_w ≳ 100 τ.

H_TT(f; T_w) for the 0.9-cycle packet at 17.48 keV (CF-B agrees to 0.3 %; steady H of Entry 65 is
0.480 at 20 GHz, 0.940 at 100 GHz; the full table is in the log):

| T_w ns | 20 GHz | 40 | 100 | 150 | H(20)/H(100) |
|---|---|---|---|---|---|
| 0.14 | 0.0063 | 0.0064 | 0.0062 | 0.0064 | 1.020 |
| 1.4 | 0.0620 | 0.0621 | 0.0621 | 0.0624 | 0.999 |
| 5 | 0.203 | 0.204 | 0.205 | 0.206 | 0.990 |
| 14 | 0.447 | 0.467 | 0.473 | 0.475 | 0.944 |
| 22 | 0.553 | 0.613 | 0.633 | 0.636 | 0.873 |
| 30 | 0.582 | 0.702 | 0.741 | 0.747 | 0.785 |
| 50 | 0.505 | 0.777 | 0.880 | 0.893 | 0.574 |
| 80 | 0.468 | 0.756 | 0.935 | 0.959 | 0.501 |
| 300 | 0.479 | 0.737 | 0.940 | 0.975 | 0.510 |

**4. The boundaries (CF-B, confirmed by TT at the same windows, 17.48 keV).**

| | H(20)/H(100) target | T_w | CF-B | TT |
|---|---|---|---|---|
| T_flat (10 % departure from flat) | 0.900 | 19.6 ns | 0.899 | 0.897 |
| T_half (halfway to steady) | 0.756 | 32.8 ns | 0.754 | 0.752 |
| T_90 (within 10 % of steady) | 0.560 | 52.2 ns | 0.560 | 0.558 |

Flat to 1 % up to 5 ns, to 6 % at 14 ns; the factor 2 is back at about 50 ns (window), i.e. a train
≥ 440 µm long at v_L. H(f; T_w) is not monotonic: at 20 GHz it overshoots to 0.58 at 30 ns and
undershoots to 0.47 at 80 ns before settling at 0.48, which is the |1 − e^{−(κ + iπ/x)T}| ripple of
period 2x Λ/v_L (75 ns at 20 GHz) damped by e^{−κT}; TT reproduces the ripple (0.501 at 80 ns against
0.511 steady). At 8 keV (1/κ = 4.8 ns) H is steady from 30 ns on and the departure from flat is 7 %.
Below 0.5 ns the ratio oscillates by ±3 % at the difference frequency: the −η image term of a real
signal truncated at T_w, of relative size 1/(2ηT_w); this is the "leakage" Entry 71 saw at 140 ps,
and it is in CF-B (TT sits on the curve), so it is a property of the estimator, not noise.

**5. What this changes.**
- Entry 71's wording "flat to ±4 % (140 ps, leakage)" is right in kind; the ±4 % is the image term
  plus the packet's group delay, both first order and both in CF-B. Entry 71's consequences stand.
- The paper's H is now H(f; T_w) = CF-B, one formula from the ps packet (H flat, estimator = kinematic
  law) to the steady train (Entry 65), with T_flat, T_half, T_90 as the design numbers for the
  long-train experiment. The packet-free short form is a stated approximation (T_w ≳ 100 τ).
- A time-resolved measurement whose window is 20–50 ns at Si 400 / 17.48 keV sits on the steep part
  of the curve: there the crystal factor depends on the window length itself, which a measurer must
  know to 10 % to keep H to 5 % (slope of the right panel of the figure).

**Checks and scope.** Symmetric Bragg, σ, normal absorption, plane wave, quasi-static, semi-infinite;
first order in S (ν ≪ 0.1); packet entering at t = 0 from the surface. Windows > 14 ns put the packet
front 10–220 Λ (0.1–2.5 mm) deep, so a real crystal must be that thick and the quasi-static
assumption (v_L/c ≪ 1) still holds; the plane-wave assumption is untouched by the depth. Anomalous
absorption (T2c) and asymmetry not included. Class: **confirmed** (CF-B) and **falsified** (CF-A as
a universal form); prior art for the windowed first-order form not searched (T2d).

---

## 2026-10-23 — Entry 76: does the running-mean degeneracy class survive larger S? (T2b, Entry 72 iii)

**Why.** Entry 72 established at S = 1.5 that order n of the transient reads only the t_tr running
mean of J_n(S A), so two histories with the same running mean are degenerate at that order (TT
difference 9e-7). Entry 72 (iii) left open whether the class breaks when S grows, where the
independent-satellite coupling is stronger and the ramp crosses Bessel extrema. BACKLOG T2b.
Script `explore/nb/e76_degeneracy_S.py`.

**Pre-registered rule (written before the run).** Laue, α = 0.5, α_ac = 40 (R = 0.0125, the
Entry 67 and 72 conditions), no absorption, 16 carrier phases at fixed envelope time, η = +q and +2q
in one call, 151 times on [0, 1.5] t_tr, `required_nz(α, α_ac, S, 2q)`. Scan S = 1.5, 2, 3, 4, 5, 6.
For each S two constructions, both from the linear ramp A_1 = clip(s, 0, 1):
- (a) degenerate at order 1: J_1(S A_2) = J_1(S A_1) + δ, δ = ε sin(2π(s − a)/(b − a)) on [a, b],
  zero outside, ε = 0.1, with [a, b] chosen inside one monotone branch of J_1(S s) with margin so
  that the inversion exists and is continuous (Newton from x0 = S A_1(s)); for S ≥ 4 a second pair
  is built on the descending branch (S s between 1.84 and 5.33), where the strain is large.
- (b) degenerate at order 2: the same with J_2 on its rising branch (S s ≤ 3.05).
The window where the whole wobble sits inside [t − 1, t] is t ∈ [b, a + 1]; there Θ_n of the
degenerate order is identical by construction (checked numerically, law difference < 1e-8).
- **Broken at order n, at this S** if max over the window of |I_n^TT(history 2) − I_n^TT(history 1)|
  exceeds 5e-4, i.e. ten times the Entry 72 convergence error 5e-5. The convergence error is
  re-measured at every S (16 → 32 phases and n_z × 2 at three times); if it exceeds 5e-5 at some S,
  the threshold there is ten times the measured value and this is stated.
- **Class holds** if no S ≤ 6 and no construction exceeds the threshold.
- Side checks, not part of the verdict: the non-degenerate order's TT difference agrees with the
  law's to 1e-3 (Entry 72 ii), and |I_TT − I_law| per S is recorded as the law's own error.
- One exploratory row after the scan, labelled as such: the (a) construction at S = 3 with a larger
  coupling ratio R = 0.1 (α = 1, α_ac = 10), to see whether R rather than S is what breaks the class.
  It is not covered by the rule above.

**Setup (as run).** As pre-registered. The interval [a, b] is chosen automatically inside the
branch so that J_n(S s) ± ε stays in the branch range (a ≥ 0.02, b ≤ 0.95); history 2 is built by
Brent inversion on a 30001-point grid and interpolated. At S = 1.5 the automatic interval
[0.134, 0.759] differs from Entry 72's [0.2, 0.7]; the order-1 result (1.2e-6 against 8.8e-7) is
the same in kind. Two runs of under two minutes each (`runA` S ≤ 3, `runB` S ≥ 4 + exploratory),
140 s in all. Reproduce: `uv run python explore/nb/e76_degeneracy_S.py runA runB fig`; data
`e76_degeneracy_S.json`, log `e76_degeneracy_S.log`, figure `figures/paper/e76_degeneracy_S.png`.

**Results.** Degenerate order: max over the window of |I_n^TT(2) − I_n^TT(1)|; "other" is the
non-degenerate order in the same window (TT, with the law's value in brackets); threshold 5e-4 at
every S (the re-measured convergence error never exceeded the Entry 72 value: 5.1e-5 at S = 1.5,
falling to 1.3e-6 at S = 6).

| S | (a) same Θ_1, J_1 rising: ΔI_1 | other ΔI_2 | (a2) same Θ_1, J_1 descending: ΔI_1 | other ΔI_2 | (b) same Θ_2: ΔI_2 | other ΔI_1 | law error |
|---|---|---|---|---|---|---|---|
| 1.5 | 1.2e-6 | 2.1e-3 (2.1e-3) | — | — | 1.2e-6 | 6.6e-3 (6.6e-3) | 6.7e-5 |
| 2 | 1.1e-6 | 2.3e-3 | — | — | 1.1e-6 | 1.4e-2 | 8.2e-5 |
| 3 | 1.3e-6 | 2.2e-3 | — | — | 3.0e-7 | 2.0e-2 | 5.0e-5 |
| 4 | 2.0e-6 | 1.9e-3 | 1.9e-6 | 9.3e-3 | 5.2e-7 | 1.3e-2 | 5.3e-5 |
| 5 | 2.6e-6 | 1.4e-3 | 2.3e-6 | 1.4e-2 | 4.0e-7 | 9.2e-3 | 6.5e-5 |
| 6 | 3.0e-6 | 9.9e-4 | 1.7e-6 | 1.0e-2 | 8.3e-7 | 6.7e-3 | 7.1e-5 |

- **Verdict under the rule: the class holds for every S ≤ 6, at both orders and on both branches
  of J_1.** The largest degenerate-order difference is 3.0e-6 (S = 6, rising branch), 170 times
  below the threshold and 20 times below the law's own absolute error. Outside the window the same
  pairs differ by 2e-3 to 3e-2, so the test has power.
- The non-degenerate order's TT difference agrees with the law's to ≤ 9e-6 at every S (Entry 72 ii
  extended to S = 6), so the law predicts which order does distinguish the two histories and by how
  much: in the (a2) pairs the order-2 difference is 1e-2 while order 1 is blind to 2e-6.
- The law's absolute error stays 5e-5 to 8e-5 over S = 1.5–6 (Entry 67's 8.7e-5 at S ≤ 5 extended
  to 6), while the convergence error falls to 1e-6; the 1e-6 degenerate-order residual is therefore
  real but is 20–50 times smaller than the law's own error, i.e. the degeneracy is more exact than
  the law that predicts it: whatever corrects the law at O(1e-4) depends on the history almost
  only through Θ_n.
- Exploratory row (not under the rule), R = 0.1 (α = 1, α_ac = 10): the degenerate-order residual
  rises to 2.4e-4 (S = 3) and 2.9e-4 (S = 6, descending branch), still under 5e-4, with the law's
  error at 2e-3. Going from R = 0.0125 to 0.1 multiplies the residual by about 170 (between R² = 64
  and R³ = 512, two points only), so the coupling ratio R, not S, is what erodes the class, as
  Entry 54 found for the independent-satellite approximation itself. At R ~ 0.2 the residual would
  reach the threshold by that scaling; not run, not claimed.

**Answer to Entry 72 (iii) / T2b.** Raising S to 6 does not break the running-mean class at
R = 0.0125; at every order the transient still reads only Θ_n to 3e-6. The order-2 reading of
Entry 72 (iii) is now quantified across S: the pairs that are degenerate at order 1 differ at
order 2 by 1e-3 (S = 1.5–6, rising branch) to 1.4e-2 (descending branch), so a two-order readout
resolves them when the order-2 satellite is measured to better than that. Class: **confirmed**
(within the model; the class is a consequence of the Θ_n law, Entry 33 claim #21, not a new
mechanism).

**Scope.** Laue, no absorption, α = 0.5, R = 0.0125 (one exploratory R = 0.1), linear ramp with one
sinusoidal zero-mean wobble of 5–10 % in J_n, orders 1 and 2, η exactly at nq. Bragg not run
(Entry 30: tanh² reading is exact at resonance, so the same degeneracy is expected). Windows where
the wobble is only partly inside [t − 1, t] are not degenerate and are not claimed to be.

---

## 2026-10-23 — Entry 77: anomalous absorption χ_ih in the Bragg solver, and its effect on H (T2c)

**Why.** Entries 59, 65, 75 use normal absorption only (η → η + iκ). The closed form H(x, g) has
one crystal parameter, g = κ/2πα = χ_i0/|χ_rh|. For Si 400, X0h (Entry 20, `e20_units_x0h.py`)
gives χ_ih/χ_i0 = 0.925 at 8 keV and 0.931 at 17.48 keV, so the two Bloch modes have very different
absorption (Borrmann) and the question is whether H needs χ_ih as a second parameter at
κ/2πα ≈ 0.01. BACKLOG T2c. Script `explore/nb/e77_chi_ih.py`; src change in `tratt.core.solve_bragg`
(new keyword `kappa_h`), tests in `tests/test_core.py`.

**Model.** For a centrosymmetric crystal χ_h = χ_hbar = χ_rh + iχ_ih, so both couplings become
one complex number: c = iπα → iπα + κ_h/2 with κ_h = 2πα χ_ih/|χ_rh| (same units as κ; κ_h/κ =
χ_ih/χ_i0). The linearisation of Entries 59/75 goes through unchanged with complex c, so the
"extended closed form" is the same formula evaluated with complex c (R0 from c R0² + iη_c R0 + c = 0,
Λ = iη_c + 2cR0). Sign: with the plus sign the semi-infinite static reflectivity is higher at η < 0,
which is the Borrmann (low-absorption) edge; in the usual convention η ∝ θ − θ_B that is the
low-angle side, where the standing-wave antinodes sit between the planes. The physical statement
is convention-free: the side with the higher static reflectivity is the low-absorption side, and
H is reported on both sides.

**Pre-registered rule (written before the run).**
1. Src: `solve_bragg(..., kappa_h=0)` reproduces the present results bit for bit at kappa_h = 0
   (existing tests unchanged, plus an explicit equality test). New tests: the thick absorbing slab
   with kappa_h > 0 reaches the root |R| < 1 of c R² + i(η + iκ) R + c = 0 with complex c to 1e-5;
   the static reflectivity at η = −πα exceeds that at +πα for kappa_h > 0 and the two are equal at
   kappa_h = 0; |R|² ≤ 1 with modulation. The full test suite runs after the change ([HEAVY],
   approved by boss, BelowNormal).
2. TT check of the extended closed form: single mode at (g, x) = (0.0100, 3 / 10 / 30 / 100) with
   κ_h = 0.93 κ, as in `e65_si400.py tt`, both signs of η. Must agree to 1 %, else the extension is
   wrong and the src change is suspect (the rule is on the physics, not on this check).
3. Physics: H(f) from the extended closed form with κ_h = 0.93 κ against κ_h = 0, Si 400 at
   17.48 keV (g = 0.0100) and 8 keV (g = 0.0458), f = 5–300 GHz, both signs of η; steady and at the
   windows 14 ns and 52 ns (CF-B of Entry 75 with complex c). **If max over 20–150 GHz of
   |H(κ_h)/H(0) − 1| exceeds 1 % at g = 0.0100 (either sign), the closed form needs χ_ih as a second
   parameter**; otherwise one parameter suffices for Si 400 and the anomalous term is a stated
   ≤ 1 % correction. The band is the Chalise et al. band (Entry 57); the number is reported
   outside it too.
4. Expectation written down beforehand (kinematic estimate): −2cR0 ≈ iπ/x + κ_h/x, so the real part
   of κ − 2cR0 changes by (κ_h/κ)/x, i.e. ~3 % at 20 GHz and ~0.7 % at 100 GHz at 17.48 keV, with
   opposite signs on the two sides; the effect on H is that times κ²/(κ² + (π/x)²). This puts the
   answer near the 1 % threshold; the rule is not adjusted for it.

**Src change (f175d12).** `solve_bragg(..., kappa, kappa_h=0.0)`: `c = i pi alpha + kappa_h / 2` for
both couplings; `kappa_h < 0` rejected. Tests added to `tests/test_core.py`: bit-identical at
`kappa_h = 0` (`np.array_equal`), complex Darwin–Prins root to 1e-5 (κ_h = κ/2 so that the
low-absorption mode still hides the back face: at κ_h → κ that mode absorbs at κ − κ_h and a
thickness of 1/κ is no longer thick, which is the Borrmann effect itself, seen first as a failing
test at κ_h = 11, κ = 12), higher static reflectivity at η = −πα than at +πα for κ_h > 0 and equal
at κ_h = 0, and the negative-value check. Full suite after the change: 58 passed in 41 s
(BelowNormal, approved HEAVY). README and DESIGN.md item 4 updated (Laue still has no χ_h'').

**Reproduce.** `uv run python explore/nb/e77_chi_ih.py run ttp ttm conv60 fig` (run 2 s; ttp
45 s; ttm 50 s; conv60 72 s); data `e77_chi_ih.json`, log `e77_chi_ih.log`, figure
`figures/paper/e77_chi_ih.png`.

**1. TT check of the extended closed form (rule item 2).** Single mode, g = 0.0100 (α = 95.4,
κ = 6), κ_h = 0.93 κ, S = 0.002 and 48 RK4 steps per cycle (Entry 65 §3; the first run of this
script copied the unconverged 24-step, S = 0.02 settings of `e65 tt` and gave 0.92 and 0.68 at
x = 30 and 100, which Entry 65 §3 already says not to quote; corrected before the numbers below,
rule unchanged).

| x | η > 0: TT / closed form | κ_h effect, TT / closed form | η < 0: TT / closed form | κ_h effect, TT / closed form |
|---|---|---|---|---|
| 3 | 1.0000 | 0.9991 / 0.9989 | 0.9982 | 0.9996 / 1.0011 |
| 10 | 1.0007 | 0.9957 / 0.9964 | 1.0028 | 1.0049 / 1.0036 |
| 30 | 1.0001 | 0.9913 / 0.9919 | 1.0015 | 1.0091 / 1.0083 |
| 60 | 0.9896 (48/cycle), **0.9968 (96/cycle)** | 0.9913 / 0.9910 | 0.9890 (48/cycle) | 1.0090 / 1.0093 |

The extended closed form holds to 0.3 % once converged (x = 60: 0.9893 → 0.9966 without κ_h and
0.9896 → 0.9968 with it, 48 → 96 steps per cycle, so the 1 % at 48 steps is the h⁴ phase error of
Entry 65 §3(a), not the extension); the κ_h-specific change of the weight agrees between TT and the
closed form to 0.03–0.15 % (η > 0) and 0.08–0.15 % (η < 0). Item 2 passes.

**2. The rule (item 3): H(κ_h)/H(0) − 1 in per cent, Si 400, κ_h/κ = 0.93.**

| f GHz | 5 | 10 | 20 | 30 | 40 | 60 | 87 | 100 | 150 | 200 | 300 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 17.48 keV, η > 0, steady | −0.25 | −0.47 | −0.78 | −0.90 | −0.91 | −0.81 | −0.65 | −0.59 | −0.42 | −0.32 | −0.22 |
| 17.48 keV, η < 0, steady | +0.25 | +0.47 | +0.79 | +0.92 | +0.94 | +0.84 | +0.67 | +0.61 | +0.43 | +0.33 | +0.23 |
| 17.48 keV, η > 0, 14 ns | −1.93 | −1.84 | −0.96 | −0.64 | −0.48 | −0.32 | −0.22 | −0.19 | −0.13 | −0.09 | −0.06 |
| 17.48 keV, η > 0, 52 ns | +0.20 | −1.49 | −1.06 | **−1.33** | −1.13 | −0.81 | −0.57 | −0.50 | −0.34 | −0.25 | −0.17 |
| 17.48 keV, η < 0, 52 ns | −0.25 | +1.68 | +1.11 | **+1.39** | +1.17 | +0.83 | +0.59 | +0.51 | +0.35 | +0.26 | +0.18 |
| 8 keV, η > 0, steady | −3.62 | −3.90 | −2.78 | −2.01 | −1.54 | −1.03 | −0.69 | −0.59 | −0.37 | −0.25 | −0.14 |
| 8 keV, η < 0, steady | +3.97 | +4.46 | +3.17 | +2.30 | +1.79 | +1.24 | +0.89 | +0.79 | +0.55 | +0.44 | +0.32 |

Maximum over the 20–150 GHz band at g = 0.0100: steady 0.91 % (η > 0) and 0.94 % (η < 0);
14 ns window 0.96 % and 0.98 %; 52 ns window 1.33 % and 1.39 %. **The rule fires, by the 52 ns
window: the closed form needs χ_ih as a second parameter** (at the steady limit and at 14 ns it
would have passed by 0.02–0.09 %, i.e. the answer sits on the threshold, as the pre-registered
estimate of item 4 said it would). At 8 keV (g = 0.046) the correction is 2.8–3.2 % in the band
and up to 4.5 % at 5–10 GHz. The pre-registered estimate (κ_h/κ)/x × κ²/(κ² + (π/x)²) gives
0.9 % at 20 GHz and 0.6 % at 100 GHz at 17.48 keV, against 0.8 % and 0.6 % computed: the mechanism
is the one written down, a 1/x correction to the effective absorption with opposite signs on the
two sides of the Bragg peak.

**3. What is new in kind: a side asymmetry.** With χ_ih the crystal factor differs between the
two detunings of the same phonon frequency, H(−η)/H(+η) = 1.016–1.019 across 20–60 GHz at
17.48 keV and 1.03–1.06 at 8 keV (up to 1.087 at 10 GHz), the Borrmann edge being the stronger one.
A satellite pair measured at ±ω therefore has a ratio that does not depend on the strain spectrum
at all (|G(ω)| cancels) and reads χ_ih/χ_i0 and x directly; without χ_ih the pair is symmetric at
every window (CF-B is even in η when c is imaginary). Not pre-registered, so a statement of what
the model predicts, not a claim; it is the natural check for the second parameter in an experiment.

**Consequences.** (i) The paper's H carries two crystal parameters, g = χ_i0/|χ_rh| and
χ_ih/χ_i0, both tabulated (X0h), so no fitting is added; the second is a ≤ 1 % steady correction
at Si 400 / 17.48 keV, 1.4 % at the 52 ns window, 3–4 % at 8 keV. (ii) Entries 59, 65, 71, 75 keep
their numbers to within those percentages. (iii) The sign convention of `kappa_h` is fixed by the
physics (higher static reflectivity on the low-absorption edge, η < 0 in the code); the detuning
side in an experiment is identified by the same static asymmetry. Class: **confirmed** (the
extension, TT to 0.3 %), **rule fired** (second parameter needed at the 1 % level), side asymmetry
**predicted, unverified**.

**Scope.** Centrosymmetric crystal (χ_h = χ_hbar), symmetric Bragg, σ, semi-infinite, first order
in S, X0h ratio 0.93 taken as the same at both energies (0.925 / 0.931). The Borrmann mode's
absorption depth is 1/(κ − κ_h) ≈ 14/κ ≈ 2.6 mm at 17.48 keV, so inside the Darwin width (x ≲ 1)
the window condition of Entry 75 changes qualitatively; not computed, outside the phonon band.
Laue χ_h'' remains excluded (DESIGN.md item 4).

---

## 2026-10-23 — Entry 78: the ±η side asymmetry of the crystal factor, tested with a packet in full TT (T2f)

**Why.** Entry 77 §3 found, in the closed form, that with χ_ih the two detunings ±η of one phonon
frequency get different crystal factors, and that their ratio A(f; T_w) = W(−η)/W(+η) cancels the
strain spectrum, so it reads χ_ih/χ_i0 directly. That was a single-mode, steady statement from the
closed form. Here it is tested with the Entry 66 packet in full TT at finite windows. BACKLOG T2f
(added by tratt with this evidence). Script `explore/nb/e78_side_asymmetry.py`.

**Pre-registered rule (written before the run).** Full semi-infinite TT (Entry 75 machinery, RK4
through the packet, analytic continuation after it) with complex coupling c = iπ + κ_h/2
(α = 1 in Λ units), Si 400 at 17.48 keV (g = 0.0100) and 8 keV (g = 0.0458), Entry 66 packet
(87 GHz, 0.9 cycles) at S/100 (peak strain 1e-5), f = 20, 30, 40, 60, 87, 100, 150 GHz at both
detunings η = ±W, windows 1.4, 14, 52, 300 ns. A_TT = W_TT(−η)/W_TT(+η); A_CF from CF-B of
Entry 75 with the same complex c.
1. κ_h = 0: A_TT = 1 within 0.5 % at every (f, T_w) (the form is even in η when c is imaginary;
   the residual is the second-order term of Entry 75 §2). If not, the TT setup is suspect and
   nothing below is claimed.
2. κ_h = 0.93 κ: if |A_TT/A_CF − 1| > 0.5 % at any (f, T_w), **the prediction is rejected** as
   stated (the closed form would be missing something the packet sees). Otherwise the asymmetry
   is confirmed in the model, with its size as computed.
3. Size, reported not judged: A − 1 in the 20–150 GHz band at each window, and dA/d(κ_h/κ) from
   the closed form (the sensitivity with which the pair reads χ_ih/χ_i0).
Convergence: S/100 → S/1000 at 14 ns (A must move by less than 0.2 %), step 64 → 128 per cycle.
Amendment before the run (boss, BACKLOG T2f wording): the cancellation of the strain spectrum in A is
checked with two packet spectra, 0.9 and 10 cycles, which must give the same A to 0.5 % at κ_h = 0.93 κ.
Run order: Entry 79 (T2e) was run before this entry, on boss's instruction.
Amendment 2 before the run (physics cross-review, PR #16): the measured weight ratio W(−η)/W(+η)
also carries |R0(−η)|/|R0(+η)|, which is not 1 with χ_ih. A is therefore defined on W/|R0| (that is
H(−η)/H(+η), what the closed form predicts), and the raw weight ratio W(−η)/W(+η) is reported next
to it with the static |R0| ratio, so that the proposal states that the pair needs the static rocking
curve as a second measurement, or an estimator built on W/|R0|.


**Run.** `uv run python explore/nb/e78_side_asymmetry.py run fig` (10 s); data `e78_side_asymmetry.json`,
log `e78_side_asymmetry.log`, figure `figures/paper/e78_side_asymmetry.png`. The first run was at
S/100 (peak strain 1e-5) and failed the pre-registered S-convergence check (A moved by 0.54 % to
S/1000) and, by the same amount, rule 1 (0.60 %): the ±η ratio is sensitive to the second-order term
at a level the single-side H of Entry 75 (0.2–0.3 %) is not. The convergence check exists for this;
the run was repeated at S/1000 (peak strain 1e-6), with S/1000 → S/10000 moving A by 0.054 % and the
step doubling by < 1e-5. Rules unchanged. Numbers below are the S/1000 run.

**Results.** A = H(−η)/H(+η) − 1 in per cent, κ_h = 0.93 κ; TT (CF-B in brackets):

| T_w | 20 GHz | 30 | 40 | 60 | 87 | 100 | 150 |
|---|---|---|---|---|---|---|---|
| 17.48 keV, 1.4 ns | 0.26 (0.22) | 0.19 (0.14) | 0.15 (0.11) | 0.09 (0.07) | 0.05 (0.05) | 0.03 (0.04) | −0.03 (0.03) |
| 17.48 keV, 14 ns | 2.00 (1.95) | 1.35 (1.30) | 1.02 (0.98) | 0.67 (0.65) | 0.45 (0.45) | 0.38 (0.39) | 0.20 (0.26) |
| 17.48 keV, 52 ns | 2.24 (2.20) | 2.80 (2.75) | 2.36 (2.32) | 1.68 (1.65) | 1.17 (1.17) | 1.01 (1.02) | 0.62 (0.68) |
| 17.48 keV, 300 ns | 1.62 (1.58) | 1.89 (1.84) | 1.91 (1.87) | 1.69 (1.67) | 1.34 (1.34) | 1.19 (1.21) | 0.80 (0.86) |
| 8 keV, 14 ns | 5.83 (5.79) | 3.95 (3.90) | 2.97 (2.93) | 1.98 (1.96) | 1.35 (1.35) | 1.16 (1.17) | 0.72 (0.78) |
| 8 keV, 300 ns | 6.17 (6.12) | 4.44 (4.39) | 3.42 (3.38) | 2.31 (2.29) | 1.59 (1.59) | 1.37 (1.38) | 0.86 (0.92) |

- Rule 1 (κ_h = 0): max |A_TT − 1| = 0.060 % at both energies (second-order residual at S/1000).
  Passed.
- Rule 2 (κ_h = 0.93 κ): max |A_TT/A_CF − 1| = 0.060 % at both energies, every window, every f.
  **The prediction is confirmed in the model.**
- Spectrum cancellation: the 0.9-cycle and 10-cycle packets give the same A to 0.105 % (both
  energies). The strain spectrum cancels, as claimed.
- The |R0| ratio: |R0(−η)|/|R0(+η)| = 1.0000 at every f here (x = 27–205, kinematic tail, where
  R0 ≈ −c/(iη_c) is even in η to O(κ_h/2η) ~ 1e-4), so the raw weight ratio equals A to 1e-4 in
  this band. The physics-review caveat applies inside and near the Darwin width (x ≲ a few), where
  the static asymmetry is 10 % (Entry 77 test, η = ±πα); a proposal that reads χ_ih from a ±ω pair
  in the phonon band needs no second measurement at Si 400, but must say so with the x range.
- Size and sensitivity: at 17.48 keV A − 1 reaches 1.9 % (steady, 30–40 GHz) and 2.8 % at the
  52 ns window (the window ripple of Entry 75 is not even in η); at 8 keV 6.2 % at 20 GHz.
  dA/d(χ_ih/χ_i0) = 0.017–0.020 (17.48 keV, 20–60 GHz) and 0.068 (8 keV, 20 GHz): a 0.1 % ratio
  measurement at 8 keV resolves χ_ih/χ_i0 to 0.015, i.e. the Debye–Waller factor of the
  reflection; at 17.48 keV to 0.05.
- At short windows the asymmetry vanishes (0.03–0.3 % at 1.4 ns), as it must: it is a property of
  the long-window crystal factor, zero when H is flat.

**Answer to T2f.** Confirmed within the model: with χ_ih the two detunings of one phonon frequency
carry different crystal factors, the ratio is independent of the strain spectrum (two packets agree
to 0.1 %) and of |R0| in the phonon band, and it is a direct measure of χ_ih/χ_i0 with sensitivity
0.02–0.07 per unit. Class: **confirmed** (full TT, two energies, four windows, two packets, no fit);
**not observed**; prior art for a ±ω Borrmann asymmetry of acoustic satellites not searched (T2d).

**Scope.** Symmetric Bragg, σ, centrosymmetric, semi-infinite, linear in S (S/1000), packet from the
surface, x ≥ 27. Inside the Darwin width the |R0| ratio must be carried (Entry 77 §3 caveat). The
asymmetry rides on the long-window H, so it needs the same window condition as Entry 75.

---

## 2026-10-23 — Entry 79: is the coupling ratio R what breaks the running-mean degeneracy class? (T2e)

**Why.** Entry 76 held the class to S = 6 at R = 0.0125 and found, in one exploratory row at
R = 0.1, a residual 170 times larger. BACKLOG T2e (boss): raise R and find where the class breaks.
Script `explore/nb/e79_degeneracy_R.py` (reuses `e76_degeneracy_S.one_case`).

**Pre-registered rule (written before the run).** Laue, no absorption, α_ac = 10 fixed (ten acoustic
periods in the slab, satellites resolved at every α used), α = 0.125, 0.5, 1, 2, 3, i.e.
R = α/α_ac = 0.0125, 0.05, 0.1, 0.2, 0.3; constructions (a) S = 3 (J_1 rising) and (a2) S = 6
(J_1 descending) of Entry 76, 16 phases, η = q and 2q, 151 times, convergence re-measured per case
(16 → 32 phases, n_z × 2).
- **R-boundary exists** if, at some R ≤ 0.3, the degenerate-order TT difference in the window
  exceeds 5e-4 (or ten times the measured convergence error if that is larger); the smallest such
  R is the boundary and the law's own error |I_TT − I_law| at that R is reported next to it.
- **R-independent** if no R ≤ 0.3 exceeds it.
- Reported, not judged: the power law of the residual against R over the five points, and whether
  the residual stays below the law's own error (Entry 76 found 20–50× below at R = 0.0125).

**Results** (`uv run python explore/nb/e79_degeneracy_R.py run fig`, 33 s; data
`e79_degeneracy_R.json`, log `e79_degeneracy_R.log`, figure `figures/paper/e79_degeneracy_R.png`).
Degenerate-order TT difference in the window (same Θ_1), threshold 5e-4 (the measured convergence
error never exceeded 3e-5):

| R (α) | (a) S = 3: ΔI_1 | law error | (a2) S = 6: ΔI_1 | law error | other order ΔI_2, TT vs law |
|---|---|---|---|---|---|
| 0.0125 (0.125) | 5.4e-6 | 7.8e-5 | 4.8e-6 | 9.0e-5 | agree to 4e-6 |
| 0.05 (0.5) | 6.1e-5 | 6.1e-4 | 7.7e-5 | 1.1e-3 | 7e-5 |
| 0.1 (1) | 2.4e-4 | 2.3e-3 | 2.9e-4 | 2.1e-3 | 3e-4 |
| 0.2 (2) | **3.3e-3** | 3.2e-2 | **1.5e-3** | 1.7e-2 | 1e-3 / 2e-3 |
| 0.3 (3) | 1.8e-2 | 1.1e-1 | 4.2e-3 | 8.0e-2 | 5e-3 / 1e-2 |

- **Verdict under the rule: an R-boundary exists, at R = 0.2** for both constructions (R = 0.1
  holds, 0.2 breaks; the boundary lies between). The law's own error at R = 0.2 is 1.7–3.2e-2,
  ten times the degeneracy residual.
- Power law of the residual: R^2.5 (a) and R^2.1 (a2) over the five points. At every R the
  residual stays below the law's absolute error by a factor 7–25: the running-mean reading remains
  more exact than the law that predicts it, all the way to R = 0.3 where the law is off by 10 %.
- Not a function of R alone: at R = 0.0125 the residual is 5.4e-6 here (α = 0.125, α_ac = 10) and
  1.3e-6 in Entry 76 (α = 0.5, α_ac = 40, same S = 3 construction). The scan fixes α_ac; a scan at
  fixed α would separate the two dependences (not run).

**Answer to T2e.** The variable that erodes the class is the coupling ratio, as the Entry 76
exploratory row suggested; the class is exact to 5e-4 for R ≤ 0.1, which is the same bound as the
independent-satellite rule of Entry 19, and fails by 3e-3 at R = 0.2 while the law is already off
by 3 %. For the paper: "order n reads the running mean of J_n(S A)" is stated with R ≤ 0.1, the
regime of every poster case. Class: **confirmed with a boundary** (within the model, Laue, two
constructions, orders 1–2).

**Scope.** Laue, no absorption, α_ac = 10, α = 0.125–3 (Pendellösung crossings at α ≥ 1 do not
change the picture), S = 3 and 6, one wobble shape.

---

## 2026-10-23 — Entry 80: corrections to Entry 73 after the physics cross-review (u0 text; harmonic column withdrawn for the shear geometry)

**Why.** Physics cross-review of PR #16 (cb86016) raised two points on Entry 73: (a) §1 says
u0 = 7.6 pm × S; (b) the second-harmonic column uses the cumulative growth law of a longitudinal
wave, while the geometry is a pure shear wave along ⟨100⟩. Both recomputed here; nothing run but
constants (`uv run python -c`, numbers below).

**(a) u0.** u0 = S d/2π with d(440) = 96.01 pm gives **15.28 pm × S**. The table of Entry 73 §2
already has 15.3 (S = 1) and 91.7 (S = 6); only the sentence in §1 was wrong (it quoted the S = 0.5
row). Verdicts unchanged, since they rest on the strain.

**(b) Harmonics of a pure shear wave along a 4-fold axis.** Accepted, with the following reasons.
- Symmetry: propagation along x = [100], polarisation u_y (h = 440 in plane). A rotation by π
  about x is a symmetry of the cubic crystal and maps u_y → −u_y, so the elastic energy is even in
  the shear strain ε_xy and contains no term cubic in it: there is no collinear, cumulative second
  harmonic of the shear wave itself. This is the standard result for shear waves in isotropic solids
  (Zarembo & Krasil'nikov, Sov. Phys. Usp. 13, 778 (1971): shear waves generate no second
  harmonic, the first self-distortion is the third harmonic; to be re-read before citing) and the
  argument above extends it to the ⟨100⟩ axis of a cubic crystal.
- The quadratic term ε_xy² ε_xx drives a longitudinal component at 2ω. It is not phase-matched:
  Δk = 2ω(1/v_T − 1/v_L) gives a coherence length π/Δk = 4.8 µm at 1 GHz and 2.4 µm at 2 GHz,
  against L = 163–296 µm, so that component is bounded by (β/4) k ε L_coh ≈ 1e-4 (S = 1) to 1.3e-3
  (S = 6, 2 GHz) of the shear strain. It is also longitudinal, u ∥ x, so h·u = 0 for the in-plane
  440 reflection: it does not enter the X-ray phase at all.
- The cubic self-distortion is of order γ (kL) ε²: (kL) ε² = 4.7e-8 (S = 1, 1 GHz) to 2.5e-5
  (S = 6, 2 GHz, 30 keV), so even with a third-order constant of order 10 the third harmonic is
  ≤ 3e-4 of the fundamental.
- The convention point: for a longitudinal wave the strain second harmonic from
  u2 = β k² u1² x/8 is ε2/ε1 = (β/4) k ε1 x, not (β/2); the (β/2) form is the fluid-acoustics one
  (p2/p1 = β ε k x/2 with ε the acoustic Mach number). Entry 73's column was therefore 2× too large
  even for the longitudinal case under the solid-state β. Moot here.

**Consequence: the ε2/ε column of Entry 73 §2 is withdrawn for the shear geometry**, and the
verdict table becomes

| S | 1 GHz | 2 GHz |
|---|---|---|
| 0.5–2 | **feasible**: ε ≤ 3.3e-5, I ≤ 25 W/cm² (pulsed) | **feasible**: ε ≤ 6.6e-5, I ≤ 100 W/cm² |
| 3–4 | **limit by drive power only**: I = 56–100 W/cm² | **limit by drive power only**: I = 230–400 W/cm² |
| 6 | **limit by drive power only**: I = 230 W/cm² | **limit by drive power only**: I = 900 W/cm²; the wave itself stays sinusoidal to ≤ 1e-3 |

Entry 73 §3's "above ~10 % the Entry 35 harmonic non-ideality applies" no longer has a case in this
geometry; Entry 35 stays relevant for longitudinal geometries (Si 400 Bragg with a longitudinal
wave, the Chalise geometry), where the (β/4) form applies. Linear elasticity, attenuation, R and
β_qs statements of Entry 73 §4 stand. Drive power is now the only limit on the elastic side for all
S ≤ 6 at 1–2 GHz; it remains unverified (DotFiles#57).

**Scope.** Symmetry argument plus order estimates, not a nonlinear elastic computation; third-order
constants of Si not looked up (γ ~ 10 assumed as an order); isotropic-solid literature statement to
be re-read. The Zarembo & Krasil'nikov reference is from memory and goes on the T2d list.

---

## 2026-10-23 — Entry 81: asymmetric Bragg geometry is an exact rescaling of the symmetric solver; the geometry factor of H (T2g)

**Why.** DESIGN.md item 2 and paper_candidates §4 list asymmetric geometry (γ_0 ≠ |γ_h|) as a
solver extension to decide on before writing, because the χ_0 gauge of the Laue solver depends on
symmetry. For Bragg the χ_0 term already enters R through η → η + iκ with
κ = μL(1/γ_0 + 1/|γ_h|)/2 (docstring of `solve_bragg`), so the question is only what the two
unequal coupling terms do. Derivation first, then a TT test of it. Script `explore/nb/e81_asym_map.py`.

**Derivation.** Two-beam Bragg, depth z, D_0' = (1/γ_0)[a_0 D_0 + b e^{+iφ} D_h],
D_h' = (1/γ_h)[(a_0 + iη') D_h + b e^{−iφ} D_0], γ_h < 0, a_0 = ikχ_0/2, b = ikCχ_h/2
(centrosymmetric, one b). For R = D_h/D_0:

    R' = −c_h e^{−iφ} − i η_c R − c_0 e^{+iφ} R²,   c_h = bL/|γ_h|,  c_0 = bL/γ_0,
    i η_c = −a_0 L (1/γ_0 + 1/|γ_h|) − iη' L/|γ_h|.

With R̃ = R √(c_0/c_h) = R √(|γ_h|/γ_0) this is the symmetric equation with one coupling
c_eff = √(c_0 c_h) = bL/√(γ_0|γ_h|), i.e. α_asym = L|χ_rh|/(λ√(γ_0|γ_h|)) (the standard
asymmetric extinction length), κ as above, and the same strain phase φ. So:
- **no solver change is needed**: `solve_bragg(alpha_asym, eta, phi, kappa, kappa_h)` and
  |R|² = |R̃|² γ_0/|γ_h| give the asymmetric case exactly (the γ_0/|γ_h| factor is the usual
  asymmetry factor of the reflectivity);
- H(x, g) keeps its form, with g → g_asym = g_sym (√b + 1/√b)/2, b = γ_0/|γ_h|, since
  κ scales with (1/γ_0 + 1/|γ_h|)/2 and 2πα with 1/√(γ_0|γ_h|); x = f/f_D with f_D = v/Λ_asym;
  the window condition of Entry 75 is in units of Λ_asym/v; χ_ih/χ_i0 is unchanged.

**Pre-registered rule (written before the run).** Semi-infinite TT with the two couplings kept
distinct (c_h ≠ c_0), Entry 66 packet at S/100, Si 400 / 17.48 keV constants, b = 4 (γ_0 = 2|γ_h|),
η at ±20, ±100 GHz, against `e66.surface_R`-style symmetric integration at c_eff with R̃ rescaled:
- if max |R_asym(t) − R̃(t) √(γ_0/|γ_h|)| over t ∈ [0, 14 ns] exceeds 1e-6 (both at the same
  RK4 grid), **the mapping is wrong** and a solver extension is needed;
- otherwise the mapping is exact and the geometry factor of H is g → g (√b + 1/√b)/2, reported
  for Si 400 at b = 1, 2, 4, 10 from the closed form (H(20)/H(100) and the window numbers T_flat,
  T_half, T_90 rescaled).

**Run.** `uv run python explore/nb/e81_asym_map.py run` (11 s); data `e81_asym_map.json`, log
`e81_asym_map.log`. The asymmetry angle χ is solved from b = sin(θ_B + χ)/sin(θ_B − χ) at
θ_B = 15.14°; Λ_asym/Λ_sym = √(γ_0|γ_h|)/sin θ_B; the same peak strain 1e-5 (S is unit-free).

**Results.**

| b = γ_0/\|γ_h\| | χ | Λ_asym µm | g_asym = g (√b + 1/√b)/2 | max \|R_asym − √b R̃\| (packet grid) | W_asym/(b W̃) − 1 at 14 ns | H_TT(b)/CF-B(g_asym) worst |
|---|---|---|---|---|---|---|
| 0.5 | −5.15° | 10.84 | 0.01061 | 7e-17 | 2e-12 | 0.30 % |
| 2 | +5.15° | 10.84 | 0.01061 | 1e-16 | 1e-12 | 0.30 % |
| 4 | +9.22° | 9.12 | 0.01251 | 0 | 0 | 0.30 % |

- **Mapping exact** (rule: < 1e-6; found machine precision), so `solve_bragg` already covers
  asymmetric Bragg: call it with α_asym = L|χ_rh|/(λ√(γ_0|γ_h|)), κ = μL(1/γ_0 + 1/|γ_h|)/2,
  κ_h likewise, and multiply |R|² by γ_0/|γ_h|. The symmetric limit is the same function, so the
  bit-identity test of T2g is satisfied by construction; no src change, no new test.
- **One geometry factor** (boss's criterion, TT vs closed form 0.3 %): H keeps its form with
  g → g (√b + 1/√b)/2 and the frequency axis in units of f_D = v/Λ_asym; the 0.30 % is the Entry 75
  second-order residual at S/100, identical for all b (as it must be for an exact mapping). b and
  1/b give the same H, since the factor is symmetric under √b ↔ 1/√b.

Closed-form consequences for Si 400 at 17.48 keV (g_sym = 0.0100), 20 against 100 GHz:

| b | Λ_asym µm | g_asym | steady H(20) / H(100) | ratio | T_flat / T_half / T_90 ns |
|---|---|---|---|---|---|
| 1 | 11.54 | 0.01001 | 0.480 / 0.939 | 0.511 | 19.6 / 32.7 / 52.3 |
| 2 (or 0.5) | 10.84 | 0.01061 | 0.479 / 0.939 | 0.510 | 17.3 / 28.8 / 46.2 |
| 4 | 9.12 | 0.01251 | 0.475 / 0.938 | 0.507 | 12.2 / 20.5 / 32.7 |
| 10 | 6.48 | 0.01740 | 0.471 / 0.937 | 0.503 | 6.1 / 10.4 / 16.5 |

- The steady crystal factor is almost geometry-independent at these x (the 20 GHz suppression
  moves from 0.480 to 0.471 at b = 10), because in the kinematic tail H ≈ κ/|κ + iπ/x| depends on
  g and x through the product g·x = κ/(2π)·f/f_D ∝ κΛ/v·f, and κΛ ∝ (1/γ_0 + 1/|γ_h|)√(γ_0|γ_h|)
  changes by only (√b + 1/√b)/2 while the frequency scale changes by 1/Λ_asym — the two nearly
  cancel for the ratio.
- The **window numbers do shrink**: T_90 goes from 52 ns (symmetric) to 33 ns at b = 4 and 16.5 ns
  at b = 10 (χ = 12.5°, grazing exit at 2.7°), because both the extinction depth Λ_asym and 1/κ in
  physical units get shorter. So an asymmetric cut is a way to bring the steady H within a shorter
  delay scan, at the price of the γ_0/|γ_h| intensity factor and a grazing beam.

**Answer to T2g.** Asymmetric Bragg needs no solver extension: it is an exact rescaling, and the
crystal factor has one geometry factor, g → g (√b + 1/√b)/2, with the window condition rescaled by
Λ_asym. Class: **confirmed** (algebra + TT to machine precision; closed form at the usual 0.3 %).
The open decision in paper_candidates §4 is closed. Laue asymmetry is a different matter (the χ_0
term does not cancel there) and is not touched.

**Scope.** σ polarisation, centrosymmetric, bulk wave along the surface normal (φ = φ(z, t)
unchanged by the cut), semi-infinite; grazing-exit geometries (b ≳ 10) bring in the specular
refraction and beam-footprint effects the plane-wave model does not have.

---

## 2026-10-23 — Entry 82: beam divergence and the crystal factor: angular convolution of H (T2k)

**Why.** `notes/referee_pass.md` objection 2 names the cheapest missing computation: what a
divergent (or polychromatic) beam does to H(f; T_w). BACKLOG T2k. Script `explore/nb/e82_divergence.py`.

**Physics written down before the run.** A ray at detuning η' = η + δ rings, after the packet, at
its own frequency Im Λ(η') ≈ η' with damping Re Λ ≈ −κ (Fig. 1a, Entry 75 §1). The Fourier
weight at the nominal ω = |η| therefore collects, from a beam of angular distribution p(δ), only
the rays within the ringing linewidth: the resonant denominator Λ(η') − iW = i(δ − π/x) − κ + …
is a complex Lorentzian in δ of half-width κ, centred at the extinction shift π/x. The beam-averaged
weight is the coherent sum W_conv(ω) = |∫ p(δ) Ĩ(η + δ; ω, T_w) dδ| (the rays' intensities add; the
oscillating parts add as phasors). Expectations, to be confirmed or refuted by TT:
(i) for σ ≫ κ + π/x the integral gives π p(π/x) × (the single-ray value × |κ − 2cR0|), so
H_conv(f) → κ √(π/2)/σ, **flat in f** and reduced by κ/σ; the steady crystal factor of Entry 65
and the window law of Entry 75 are then unobservable; (ii) with a finite window the ray response
has width ~1/T_w instead of κ, so the window law survives only while 1/T_w ≳ σ, i.e. T_w ≲ 1/σ
(0.7 ns at σ = 0.3 Darwin half-widths): the effective window is min(T_w, 1/σ, 1/κ); (iii) in
angle, κ is 0.01 (17.48 keV) and 0.046 (8 keV) Darwin half-widths: 0.03 µrad and 0.4 µrad; in
energy, 2 and 4 meV. If (i) holds, the long-window H needs a beam collimated and monochromatised
to that level, or an angle-resolving detector, and the ps-packet proposal must say so.

**Pre-registered rule (BACKLOG T2k, written before the run).** Full semi-infinite TT (Entry 75
machinery), Entry 66 packet at S/1000, Si 400 / 17.48 keV, f = 20, 40, 100 GHz, windows 1.4,
14, 52, 300 ns, Gaussian angular distribution of width σ = 0.3, 1, 3 Darwin half-widths (σ_η =
2πσ in η units), on a non-uniform δ grid that resolves κ (spacing ≤ κ/3 near the centre) and
reaches ±6σ. TT angle average: W_TT,conv(ω) = |Σ_k p_k Ĩ_TT(η + δ_k; ω)|, each ray a full TT run
(vectorised over η). Closed form: the same sum with CF-B of Entry 75 evaluated at detuning η + δ_k
and Fourier frequency ω = |η|.
- If |W_TT,conv / W_CF,conv − 1| > 1 % at any (f, σ, T_w), "a divergence term is needed";
  otherwise "closed form retained" (the convolution of CF-B is the closed form with divergence).
- Reported either way: H_conv(f; T_w, σ)/H_conv(100 GHz) against T_w for each σ, with T_flat and
  T_90 where they still exist, and H_conv/H_single at 300 ns against the expectation κ√(π/2)/σ.
- Convergence: δ-grid spacing halved at one (f, σ, T_w) must change W_conv by < 0.2 %.

**Run.** `uv run python explore/nb/e82_divergence.py run fig` (13 s); data `e82_divergence.json`, log
`e82_divergence.log`, figure `figures/paper/e82_divergence.png`. The first run used a δ grid with
du = 0.1 and failed the pre-registered 0.2 % grid-convergence check by 2 %; the grid was refined to
du = 0.025 (527–711 rays per σ), where halving again changes W_conv by 6e-7 relative. Rule unchanged.

**Results.** H_conv = beam-averaged Fourier weight at ω = |η| over the kinematic law of the nominal
ray (so σ = 0 is the Entry 75 H); Si 400, 17.48 keV, S/1000.

| σ (Darwin half-widths) | T_w | H_conv at 20 / 40 / 100 GHz | H(20)/H(100) | TT / closed form |
|---|---|---|---|---|
| 0 | 1.4 / 14 / 52 / 300 ns | 0.062 / 0.448 / 0.497 / 0.480 (20 GHz) | 1.001 / 0.946 / 0.561 / 0.511 | 0.9998–1.0001 |
| 0.3 | 1.4 / 14 / 52 / 300 ns | 0.0386 / 0.0406 / 0.0406 / 0.0406 (all f equal to 0.3 %) | 0.999 / 0.998 / 0.998 / 0.998 | 0.9998–1.0001 |
| 1 | all | 0.0124 (all f, all T_w) | 1.000–1.004 | 0.9998–1.0001 |
| 3 | all | 0.0042 (all f, all T_w) | 1.007–1.013 | 0.9998–1.0001 |

- **Rule: closed form retained.** Worst |W_TT,conv/W_CF,conv − 1| = 0.021 % over all (f, σ, T_w):
  the convolution of CF-B over the beam is the closed form with divergence; no extra term.
- **The expectation (i) holds.** For σ ≥ 0.3 Darwin half-widths H_conv is flat in f to 0.3 %
  (σ = 0.3) and 1.3 % (σ = 3) at every window, and its value is κ√(π/2)/σ_η: 0.0418 / 0.0125 /
  0.0042 expected against 0.0407 / 0.0124 / 0.0042 found (3 % low at σ = 0.3, where σ_η = 1.9 is
  only 14× the Lorentzian half-width κ + π/x = 0.13 at 20 GHz).
- **Expectation (ii) holds.** The window law disappears: H_conv is the same at 1.4 ns and at
  300 ns (0.0386 vs 0.0406 at σ = 0.3, the 5 % difference being the 1/σ_η = 0.53 Λ/v = 0.7 ns
  coherence time against the 1.4 ns window). The effective window is min(T_w, 1/σ_η, 1/κ); for
  σ_η ≫ κ it is 1/σ_η, independent of f, which is why the crystal factor's shape is gone.
- **Scale (iii).** κ in Darwin half-widths is g: at 17.48 keV the ringing half-width is 0.032 µrad
  (Darwin half-width 3.18 µrad) or 2.1 meV; at 8 keV 0.38 µrad (of 8.33 µrad) or 4.4 meV. T_flat
  and T_90 exist only for σ ≲ 0.01–0.05 Darwin half-widths; at σ ≥ 0.3 the ratio never departs
  from 1 by 10 % and the boundaries do not exist.

**What this means (the honest consequence for the paper).**
1. The steady crystal factor of Entry 65, the window law and design numbers of Entry 75, the χ_ih
   correction of Entry 77 and the ±η asymmetry of Entry 78 are properties of a *single ray*: they
   are observable only with an incident beam whose spread in detuning is below the ringing width
   κ + π/x, i.e. collimated to ≲ 0.03–0.4 µrad **and** monochromatic to ≲ 2–4 meV at Si 400, or
   with an angle-resolving detector of that resolution (which removes the angular part only; the
   energy part needs the monochromator). Ordinary beams (µrad, 1e-4 bandwidth = 35 µrad-equivalent
   at a Si 111 monochromator) are 100–1000× too broad.
2. For any such ordinary beam the measured Fourier weight is the **kinematic law times a flat
   factor** κ√(π/2)/σ_η: the ringing of different rays dephases after 1/σ_η (0.7 ns at
   σ = 0.3 half-widths, 0.07 ns at 3), so the estimator "divide by |R0|/ω" is correct as it stands,
   with no crystal factor, as soon as the window exceeds the coherence time. This is the practical
   content for the ps-ultrasonics community: Chalise et al.'s flat crystal factor (Entry 71) is
   not only because their window is short but because any real beam makes it flat.
3. The flat factor κ√(π/2)/σ_η is itself a measurable: it reads the beam's spread in detuning
   (angle plus bandwidth) in units of the absorption, i.e. the coherence time of the ringing.
   Not pre-registered; stated, not claimed.
4. P1's experimental weight falls: the long-window measurement is an inelastic-scattering-grade
   (meV) experiment, not a routine delay scan. The ±η asymmetry (Entry 78) is odd in δ and averages
   to zero over a symmetric beam, so it also needs the single-ray resolution. The formation law
   (P2) is a different observable (satellite centre intensity at η = q with S ~ 1) and is not
   addressed here; its divergence sensitivity is the satellite's own rocking width (α|J_n| Darwin
   half-widths), a separate question.
5. Entries 65, 71, 75, 77, 78, 81 stand as single-ray results; what changes is their experimental
   reach, which `notes/paper_candidates.md` and `notes/experiment_proposal.md` now state.

**Class.** Confirmed (full TT beam average vs closed form 0.02 %; the flat limit from the first-
order lineshape). Prior art: that a crystal's impulse response rings at the detuning and dephases
over the beam's angular acceptance is textbook dynamical diffraction; its consequence for the
Fourier-weight estimator was not found stated in the papers read (T2d unread list applies).

**Scope.** Gaussian spread in detuning (angle and bandwidth combined into one σ), first order,
Si 400 / 17.48 keV, f = 20–100 GHz, symmetric Bragg, σ polarisation. The energy-to-detuning
conversion uses tan θ_B (symmetric Bragg). Not run at 8 keV (same formula, g = 0.046).

---

## 2026-10-23 — Entry 83: the beam condition at GHz-train frequencies, x ≈ 1–10 (T2m)

**Why.** Entry 82 was run at x = f/f_D = 27–137 (the 20–100 GHz band), where the ringing
half-width κ + π/x is 0.01–0.02 Darwin half-widths. For a GHz train x ≈ 1–5 and the extinction
term π/x dominates: in Darwin half-widths the width is g + 1/(2x) = 0.15–0.5, so an ordinary
high-resolution beam may already resolve it. BACKLOG T2m. Script `explore/nb/e83_lowx.py`
(reuses the Entry 82 functions).

**Pre-registered rule (BACKLOG T2m, written before the run).** Si 400 at 17.48 keV (g = 0.0100)
and 8 keV (g = 0.0458), symmetric Bragg, Entry 66 packet at S/1000 (its spectrum at 0.5–8 GHz is
finite and cancels in the ratio), x = 1, 2, 5, 10 (f = x f_D: 0.73–7.3 GHz at 17.48 keV, 0.52–5.2
GHz at 8 keV), windows 52 and 300 ns, Gaussian beam spread σ = 0.03, 0.1, 0.2, 0.3, 0.5, 1, 2
Darwin half-widths. For each (E, x, T_w) the table gives H_conv/H_single and the largest σ at
which |H_conv/H_single − 1| ≤ 5 % (linear interpolation between scanned σ), in Darwin half-widths,
µrad and meV.
- If that σ_max at x = 1–5 exceeds ordinary beamline values (a few µrad in angle, a few meV to
  tens of meV in energy with a high-resolution monochromator; ID28/HERIX-class), **"the GHz train
  is measurable with an ordinary beam"** and the E scores of P2 and P1+P2 are re-evaluated;
  otherwise "single-ray condition stands".
- Expectation written beforehand: a Gaussian average of a Lorentzian of half-width w leaves the
  centre value within 5 % for σ ≲ 0.25 w, so σ_max ≈ 0.25 (g + 1/(2x)) half-widths: 0.13 (x = 1),
  0.07 (x = 2), 0.03 (x = 5) at 17.48 keV, i.e. 0.4, 0.2, 0.1 µrad; at 8 keV 1.2, 0.6, 0.4 µrad and
  the corresponding meV. If this holds the relaxation is real but to ~µrad, not to "ordinary".
- Convergence: δ-grid halving at one (E, x, σ, T_w) must change W_conv by < 0.2 %; TT vs closed
  form convolution must agree to 1 % (else the Entry 82 setup fails at low x and nothing is claimed).

**Run.** `uv run python explore/nb/e83_lowx.py run fig` (17 s); data `e83_lowx.json`, log
`e83_lowx.log`, figure `figures/paper/e83_lowx.png`. TT beam average vs closed-form convolution:
worst 0.24 % (rule: 1 %); δ-grid halving changes W_conv by < 1e-7. Windows 52 and 300 ns give the
same numbers at 8 keV and differ by ≤ 2 % at 17.48 keV; the 300 ns (steady) values are quoted.

**Results.** H_conv/H_single at 300 ns (σ in Darwin half-widths):

| E | x (f) | σ = 0.03 | 0.1 | 0.2 | 0.3 | 0.5 | 1 | 2 | σ_max (5 %) |
|---|---|---|---|---|---|---|---|---|---|
| 17.48 keV | 1 (0.73 GHz) | 1.14 | 1.45 | 2.01 | 2.00 | 1.56 | 0.85 | 0.23 | < 0.03 hw (< 0.1 µrad, < 6 meV) |
| | 2 (1.46 GHz) | 1.01 | 1.21 | 1.08 | 0.83 | 0.55 | 0.29 | 0.13 | 0.043 hw = 0.14 µrad = 9 meV |
| | 5 (3.65 GHz) | 1.12 | 0.97 | 0.57 | 0.39 | 0.24 | 0.13 | 0.07 | < 0.03 hw |
| | 10 (7.3 GHz) | 1.14 | 0.56 | 0.30 | 0.21 | 0.13 | 0.06 | 0.03 | < 0.03 hw |
| 8 keV | 1 (0.74 GHz) | 1.01 | 1.07 | 1.28 | 1.24 | 0.96 | 0.51 | 0.14 | 0.076 hw = 0.63 µrad = 7 meV |
| | 2 (1.49 GHz) | 1.01 | 1.13 | 0.99 | 0.78 | 0.54 | 0.29 | 0.12 | 0.053 hw = 0.44 µrad = 5 meV |
| | 5 (3.72 GHz) | 1.04 | 0.86 | 0.55 | 0.40 | 0.25 | 0.14 | 0.08 | 0.037 hw = 0.31 µrad = 3.6 meV |
| | 10 (7.44 GHz) | 0.95 | 0.59 | 0.35 | 0.25 | 0.16 | 0.08 | 0.04 | 0.031 hw = 0.26 µrad = 2.9 meV |

(hw = Darwin half-width: 3.18 µrad / 206 meV at 17.48 keV, 8.33 µrad / 96 meV at 8 keV. "< 0.03"
means the smallest scanned σ already departs by more than 5 %; the scan does not resolve below it.)

- **Verdict under the rule: the single-ray condition stands.** The relaxation at low x is real
  but modest: at 8 keV the 5 % tolerance allows 0.6 µrad / 7 meV at x = 1 and 0.4 µrad / 5 meV at
  x = 2, against 0.38 µrad / 4.4 meV at x ≥ 27 (Entry 82); at 17.48 keV it is ≤ 0.14 µrad / 9 meV.
  None of this is "a few µrad"; in energy, 5–9 meV is a high-resolution monochromator, not an
  ordinary one. The E scores of P2 and P1+P2 are therefore not raised on this account.
- The pre-run expectation σ_max ≈ 0.25 (g + 1/(2x)) was too optimistic by about 2× (0.136 vs
  0.076 at x = 1, 8 keV) and wrong in form at x ≲ 2: there the beam average first *raises* the
  weight (up to 2.0× at 17.48 keV, 1.3× at 8 keV, for σ = 0.2–0.3 hw) before suppressing it,
  because the single-ray resonance is centred at the extinction shift π/x away from the nominal
  detuning (0.5 Darwin half-widths at x = 1): a beam of that width lands on the peak that the
  nominal ray misses. The 5 % criterion is therefore violated by overshoot, not only by
  dephasing, and the beam-averaged factor at x ~ 1–2 is a strong, non-monotonic function of σ.
- At σ ≥ 1 hw the flat law of Entry 82 is approached only for x ≳ 10 (17.48 keV: H_conv = 0.0069
  at x = 10 against κ√(π/2)/σ_η = 0.0063 at σ = 2 hw); at x = 1–5 the beam-averaged factor is
  smaller still (0.0009–0.0073), so the "ordinary beam → kinematic law × flat factor" statement
  holds for the ps-packet band and not for a GHz train, where the factor depends on x and σ
  together.

**What this changes.** (i) For the GHz-train (P2) apparatus the beam condition is 0.3–0.6 µrad
and 4–7 meV at 8 keV, somewhat easier than for the ps band but still a high-resolution beam;
`notes/experiment_proposal.md` §0 is corrected (the "ordinary high-resolution beam" sentence is
withdrawn). (ii) DotFiles#56 (poster Q&A on the withdrawn Si 400 factor 2) gets the stronger
reason: with an ordinary beam there is no crystal factor at all in the ps band, and at GHz the
beam-averaged factor is not the single-ray one either. (iii) The satellite-centre observable of
P2 (formation law at η = q, S ~ 1) was not computed here; its divergence sensitivity is a separate
question (the satellite's own rocking width, α|J_n| Darwin half-widths) and stays open.

**Class.** Confirmed (full TT beam average, two energies, four x, seven σ, two windows, no fit).
**Scope.** As Entry 82; x = 1 sits at the Darwin edge where the kinematic normalisation is only a
common factor; the Entry 66 packet's spectrum at 0.7–7 GHz cancels in the ratio; symmetric Bragg
Si 400, while the poster's train is Si 440 Laue (not computed; the Laue satellite is a different
observable).

**Addendum to Entries 82–83 after the physics cross-review (PR #16, same day).** (a) Independent
reproduction of Entry 82 agreed (0.0406 / 0.0124 / 0.561). (b) The phasor sum is exact for any
transverse coherence (reflection diagonal in k_x, Parseval over exit angles), stronger than the
"intensities add" argument used above; it fails only for a partly angle-resolving detector with a
coherent beam, or a laterally non-uniform front. (c) **Correction:** the quadrature sum of angle and
bandwidth ("1e-4 bandwidth = 35 µrad-equivalent") assumes them uncorrelated; after a crystal
monochromator they are DuMond-correlated and only the projection onto η counts: dispersive
(Si 111 DCM + Si 400 sample) is as written, a non-dispersive (+,−) matched arrangement cancels the
bandwidth in η and leaves the first crystal's intrinsic width, and an asymmetric (+,−) collimator
with |b| ≪ 1 gives σ_η ≈ √|b| half-widths with routine optics. σ_η must be computed from the optics.
The single-ray numbers (0.03–0.4 µrad, 2–4 meV) stand as the requirement on σ_η. (d) The 2 meV
requirement is transform-compatible with ≥ 0.9 ps pulses; synchrotron bunches of 50–100 ps cannot
resolve 20 GHz regardless. (e) **Same phenomenon, two normalisations:** physics' exploratory
CF-B value "0.756 at x = 2, σ = 0.3, g = 0.010" is H_conv(x = 2)/H_conv(flat) = p(π/x)/p(0) ≈
exp[−(π/x)²/2σ_η²] = 0.71, and Entry 83 gives exactly that: 0.830 × 0.037 / 0.0406 = 0.756. Both
say: an ordinary beam does not see the single-ray H, but in the GHz regime (x ≲ 10, f ≲ f_D/(2σ))
it sees a frequency-dependent, non-monotonic beam-averaged factor that carries the beam's
η-distribution at the extinction shift π/x; this factor is computable (CF-B convolution, TT to
0.24 %) and falsifiable once σ_η is known from the optics. It is a pre-registration candidate
(not run as such; the Entry 83 scan is its first map), and it means P1 + P2 still share an
observable in the GHz regime even though the single-ray H needs an IXS-grade beam.

---

## 2026-10-23 — Entry 84: beam divergence and the formation law at the satellite centre (T2n)

**Why.** Entries 82–83 settled the beam condition for the Fourier-weight estimator (P1). The
poster's observable (P2) is different: the centre intensity of the n-th satellite, I_n(t̃) =
sin²(πα|J_n(S)|t̃) in Laue, tanh² in Bragg (Entries 1, 9; Fig. 4b), a step-front drive at S ~ 1.
Its angular scale is the satellite's own rocking width, w_n = 2πα|J_n(S)| in η (the Darwin-like
half-width with coupling α_n = α|J_n|), which is |J_n(S)| = 0.3–0.6 main Darwin half-widths, i.e.
0.5–1 µrad at Si 440 / 17.48 keV. BACKLOG T2n. Script `explore/nb/e84_satellite_divergence.py`.

**Pre-registered rule (BACKLOG T2n, written before the run).** The four poster cases
(S, α, n) = (0.5, 2, 1), (1.5, 2, 1), (2, 4, 2), (1.5, 4, 1), α_ac = 20, step front, Laue and Bragg,
t̃ = 0–1 (61 points + the poster's offset 0.0037), full TT `rocking_cube` on an η grid
nq + δ with spacing w_n/40 to ±5σ_η, Gaussian weights; σ_η = 0.1, 0.3, 1 × w_n. Metrics on the
first maximum of I_n(t̃) (parabolic interpolation): its time (the "period") and its height (the
"amplitude"), beam-averaged against single-ray.
- **"P2 is measurable with an ordinary beam"** if at σ_η = 1 × w_n both |Δperiod| ≤ 1 % and
  |Δamplitude| ≤ 5 % in every case and geometry; otherwise the table of the largest σ_η meeting
  both, with w_n in µrad for Si 440 / 17.48 keV (Laue) and Si 400 / 8 keV (Bragg).
- Check, not verdict: the beam average of the independent-satellite closed form (two-beam slab of
  thickness t̃ with coupling α_n at detuning δ, computed with the lossless solvers at φ = 0 by the
  scaling α → α_n t̃, η → δ t̃) must agree with the TT beam average to 2 %, else the independent-
  satellite picture fails off-centre and the closed form is not used.
- Expectation written beforehand: the satellite rocking curve at the first maximum has half-width
  ≈ w_n, so σ_η = 0.3 w_n should cost ~5–10 % in amplitude and shift the maximum earlier by a few
  per cent, and σ_η = w_n should halve the amplitude: the rule is expected to fail at 1 × w_n and
  the condition to land at σ_η ≈ 0.2–0.3 w_n, i.e. 0.1–0.3 µrad — the same sub-µrad class as P1.
- Reported, not pre-registered: for σ_η ≫ w_n the beam average becomes the integrated satellite
  population, P_n(t̃) ∝ ∫_0^{2πα_n t̃} J_0(x) dx in Laue (the classic integrated Laue intensity of a
  slab of coupling α_n t̃), which is a divergence-free observable with its own formation law; its
  TT beam average at σ_η = 3 w_n is compared with that integral.

**Run.** `uv run python explore/nb/e84_satellite_divergence.py run fig` (131 s, one run; eight
case–geometry pairs × five σ); data `e84_satellite_divergence.json`, log
`e84_satellite_divergence.log`, figure `figures/paper/e84_satellite_divergence.png`. Single-phase
samples on the poster's offset grid (no carrier-phase average), as in Fig. 4b.

**Results.** First maximum of the beam-averaged centre intensity against the single ray:
Δperiod = time shift, Δamp = height change. w_n = 2πα|J_n(S)|; |J_n(S)| is w_n in main Darwin
half-widths; the last two columns give w_n in angle.

| geometry | (S, α, n) | σ_η = 0.1 w_n: Δperiod / Δamp | 0.3 w_n | 1 w_n | 3 w_n | largest σ_η meeting both (scan) | w_n, Si 440 / 17.48 keV | w_n, Si 400 / 8 keV |
|---|---|---|---|---|---|---|---|---|
| Laue | (0.5, 2, 1) | +0.1 % / −0.3 % | +0.7 % / −2.6 % | −7 % / −25 % | −7 % / −62 % | 0.3 w_n | 0.43 µrad | 2.0 µrad |
| Laue | (1.5, 2, 1) | +0.8 % / −0.6 % | +1.1 % / −6.7 % | −14 % / −40 % | −19 % / −73 % | 0.1 w_n | 0.99 µrad | 4.7 µrad |
| Laue | (2, 4, 2) | −1.1 % / −1.3 % | −5.9 % / −7.7 % | −19 % / −38 % | −32 % / −72 % | < 0.1 w_n | 0.62 µrad | 2.9 µrad |
| Laue | (1.5, 4, 1) | +11 % / −0.8 % | −2.9 % / −7.8 % | −12 % / −39 % | −11 % / −71 % | < 0.1 w_n (metric, see below) | 0.99 µrad | 4.7 µrad |
| Bragg | (0.5, 2, 1) | 0.0 % / −0.2 % | −0.2 % / −1.4 % | −1.1 % / −13 % | −1.3 % / −52 % | 0.3 w_n | 0.43 µrad | 2.0 µrad |
| Bragg | (1.5, 2, 1) | −0.1 % / −0.1 % | −0.9 % / −1.0 % | −1.0 % / −20 % | −1.0 % / −62 % | 0.3 w_n | 0.99 µrad | 4.7 µrad |
| Bragg | (2, 4, 2) | 0.0 % / −0.1 % | −0.2 % / −0.7 % | −0.2 % / −21 % | +0.7 % / −63 % | 0.3 w_n | 0.62 µrad | 2.9 µrad |
| Bragg | (1.5, 4, 1) | 0.0 % / 0.0 % | −0.2 % / −0.1 % | −2.4 % / −21 % | −9 % / −60 % | 0.3 w_n | 0.99 µrad | 4.7 µrad |

- **Verdict under the rule: the condition table applies.** At σ_η = 1 w_n the amplitude falls by
  13–21 % (Bragg) and 25–40 % (Laue) in every case, far outside 5 %. The largest scanned σ_η that
  meets both criteria is 0.3 w_n in Bragg (amplitude −0.1 to −1.4 %, period ≤ 0.9 %) and
  0.1–0.3 w_n in Laue (amplitude −0.3 to −1.3 % at 0.1 w_n, −2.6 to −7.8 % at 0.3 w_n). In
  angle that is **0.1–0.3 µrad at Si 440 / 17.48 keV (the poster's Laue case) and 0.6–1.4 µrad at
  Si 400 / 8 keV Bragg**: the same sub-µrad class as the single-ray crystal factor of Entries 82–83,
  as the pre-run expectation said (0.2–0.3 w_n).
- The "period" metric (time of the first maximum, parabolic interpolation on a 1/60 grid) is noisy
  where the maximum is broad: the +11 % at σ_η = 0.1 w_n for Laue (1.5, 4, 1) is a detection
  artefact (amplitude −0.8 % there). A post-hoc metric computed from the saved traces, the
  half-rise time t_1/2, gives ≤ 0.5 % at 0.1–0.3 w_n in Laue (−2.9 % for (1.5, 4, 1) at 0.3),
  +1.6–4.6 % in Bragg at 0.3 w_n (the tanh² rise steepens off-centre), −3 to −29 % at 1 w_n, and
  does not change the verdict. Labelled post hoc; the rule's metric stands in the table.
- Bragg is the more tolerant geometry: the tanh² saturation makes the plateau height insensitive
  to the beam (−1 % at 0.3 w_n) while Laue's sin² maximum is a narrow feature in η.
- Check (closed form): the beam average of the independent-satellite slab form differs from the
  TT beam average by 3–15 % in absolute intensity already at σ_η = 0, i.e. the single-phase
  samples carry the fringe contamination the poster avoided by its grid offset (Entry 2); the
  2 % tolerance fails and, as pre-registered, the closed form is not used for the verdict (TT only).
- Reported: at σ_η = 3 w_n the Laue beam average is 12–15 % away from the ∫_0^{2πα_n t̃} J_0 law
  (shape, after scaling); the integrated-population observable is divergence-tolerant but its
  formation law is not that integral at this σ (the Gaussian is not flat over the satellite);
  not pursued.

**Answer to T2n.** The poster's centre-intensity formation law needs σ_η ≲ 0.3 w_n, i.e.
0.1–0.3 µrad in the poster's Si 440 / 17.48 keV Laue geometry (0.6–1.4 µrad at Si 400 / 8 keV
Bragg), so **P2 is not measurable with an ordinary beam as a centre observable** either; its
beam-tolerant form is the integrated satellite population, whose formation law is a separate,
uncomputed item. The E scores of P2 and P1 + P2 stay at 3 with this condition added; the shared
GHz observable of Entries 82–83 (the beam-averaged factor) remains the ordinary-beam route.
Class: **confirmed** (full TT beam average, eight case–geometry pairs, no fit). Prior art: the
sensitivity of a Pendellösung-type oscillation to beam divergence is textbook; not searched for
satellites.

**Scope.** Step front, α_ac = 20, single carrier phase on the offset grid, Gaussian σ_η, t̃ ≤ 1,
Laue lossless and Bragg lossless (no κ), first maximum only; the integrated population and the
ramp regime (Entry 67) not covered.

---

## 2026-10-23 — Entry 85: one GHz-train case in one table, Si 440 at 17.48 keV, 1–2 GHz (T2o)

**Why.** BACKLOG T2o: move the poster's conditions (Entry 62: symmetric Laue, L = 3Λ_P, bulk shear
wave) to Si 440 / 17.48 keV / 1–2 GHz and put S, ν, the window numbers and the beam conditions in
one table from the existing formulas and constants. No new physics, nothing pre-registered.
Script `explore/nb/e85_si440_case.py` (1 s), data `e85_si440_case.json`, log `e85_si440_case.log`.
Constants: X0h (Entry 20): Λ_P = 54.34 µm, θ_B = 21.68°, |χ_rh| = 1.213e-6, χ_i0 = 1.606e-8
(g = 0.0132), d = 96.0 pm; v_T = 5846 m/s, c44 = 79.6 GPa (Entry 73). L = 163 µm, t_tr = 27.9 ns;
Darwin half-width 1.77 µrad = 78 meV; symmetric-Bragg extinction length of the same reflection
21.6 µm, f_D(Bragg, shear) = 0.271 GHz, so 1–2 GHz is x = 3.7–7.4.

**Table (1 GHz; 2 GHz differs only where marked).** u0 = S d/2π, ε = S d/Λ_ac (Entry 73),
I = v_T c44 ε²/2 (time-averaged, travelling wave), α_ac = L/Λ_ac, R = Λ_ac/Λ_P, ν = πα|J_1(S)|
(first-order satellite amplitude of the filled Laue crystal, Entry 65's form with ℓ = L), τ_1 =
t_tr/(2α|J_1(S)|) (first maximum of the formation law), β_qs = 2πfL/c (Entry 62); beam columns:
the centre-law condition 0.3 w_1 = 0.3|J_1(S)| Darwin half-widths (Entry 84) and the single-ray
ringing width g + 1/(2x) half-widths (Entries 82–83).

| S | u0 pm | ε (1 / 2 GHz) | I W/cm² (1 / 2 GHz) | ν | τ_1 / t_tr | τ_1 ns | beam, centre law µrad | beam, ringing µrad (1 / 2 GHz) |
|---|---|---|---|---|---|---|---|---|
| 0.01 | 0.15 | 1.6e-7 / 3.3e-7 | 0.002 / 0.006 | 0.047 | 33 | 930 | 0.004 | 0.26 / 0.14 |
| 0.02 | 0.31 | 3.3e-7 / 6.6e-7 | 0.003 / 0.01 | 0.094 | 17 | 465 | 0.01 | |
| 0.1 | 1.5 | 1.6e-6 / 3.3e-6 | 0.06 / 0.25 | 0.47 | 3.3 | 93 | 0.03 | |
| 0.3 | 4.6 | 4.9e-6 / 9.9e-6 | 0.56 / 2.3 | 1.40 | 1.12 | 31 | 0.08 | |
| 0.5 | 7.6 | 8.2e-6 / 1.6e-5 | 1.6 / 6.3 | 2.28 | 0.69 | 19 | 0.13 | |
| 1 | 15.3 | 1.6e-5 / 3.3e-5 | 6.3 / 25 | 4.15 | 0.38 | 10.6 | 0.23 | |
| 1.5 | 22.9 | 2.5e-5 / 4.9e-5 | 14 / 56 | 5.26 | 0.30 | 8.3 | 0.30 | |
| 2 | 30.6 | 3.3e-5 / 6.6e-5 | 25 / 100 | 5.44 | 0.29 | 8.1 | 0.31 | |

α_ac = 27.9 / 41.8 / 55.8 and R = 0.108 / 0.072 / 0.054 at 1 / 1.5 / 2 GHz; β_qs = 3.4e-3 to
6.8e-3 (Entry 62). Window law for the same reflection in symmetric Bragg [corrected, see the
addendum: the shear wave along the normal has h·u = 0 in that geometry, so the Bragg counterpart
uses a longitudinal [110] wave, v_L = 9136 m/s, x = 2.36 / 4.73 at 1 / 2 GHz] (CF-B, Entry 75,
ratio H(1 GHz)/H(2 GHz), steady 0.484 with H = 0.059 / 0.123): T_flat 2 ns, T_half 5 ns,
T_90 7 ns; 1/κ = 260 µm = 28 ns at v_L. (First version, withdrawn: shear v_T, x = 3.7 / 7.4,
T_flat 7 / T_half 12 / T_90 17 ns.)

**Reading.**
- The two regimes of the P1 + P2 experiment sit at very different S: the linear regime (ν ≤ 0.1,
  Entry 68) of a *filled* Laue crystal at α = 3 needs S ≤ 0.02, i.e. strain ≤ 7e-7 and
  3–10 mW/cm², with a formation time of 0.5–1 µs that the transit cuts off (the satellite never
  forms; it stays kinematic, which is the point); the poster's formation regime is S = 0.5–2,
  strain 1–7e-5, 1.6–100 W/cm², τ_1 = 8–19 ns inside t_tr = 28 ns. Both are reachable with the
  same drive at two amplitudes, which is the "one apparatus" of the recommendation.
- The window numbers of the Bragg counterpart (longitudinal wave: 2–7 ns) are well inside one
  transit, so a long-train measurement at Si 440 reaches the steady crystal factor within t_tr;
  but the beam condition for that single-ray factor is 0.21–0.39 µrad (ringing width
  g + 1/(2x) at x = 2.4–4.7), ≈ 0.1–0.2 µrad at the 5 % tolerance by the Entry 83 scaling, and
  for the centre formation law 0.08–0.31 µrad (0.3 w_1): the GHz train at Si 440 is a sub-µrad
  measurement on both counts, as Entries 83–84 found for Si 400.
- The ordinary-beam observables for this case are therefore the beam-averaged factor of Entries
  82–83 (x = 3.7–7.4, non-monotonic in σ_η; computable once σ_η is known) and the integrated
  satellite population (T2p).

**Scope.** Table only; no run beyond the 1 s closed-form window numbers; Laue absorption is not in
the solver (the Laue contrast 1/cosh(2πα''J_1) of Entry 20 is 0.99 here); transducer output
unverified (DotFiles#57).

---

## 2026-10-23 — Entry 86: does the integrated satellite population lift the beam condition? (T2p)

**Why.** Entry 84: the centre formation law needs σ_η ≲ 0.3 w_n. A detector that integrates the
whole satellite channel measures the population instead, which for a beam much wider than the
satellite is the rocking-curve integral and should not care about σ_η. BACKLOG T2p. Script
`explore/nb/e86_population_divergence.py`.

**Observable and detector.** The n-th channel is |η − nq| ≤ q/2 (q = 2πα_ac = 126 at α_ac = 20,
i.e. 6.7 main Darwin half-widths, 12 µrad at Si 440 / 17.48 keV); the detector accepts the whole
channel and nothing else (an angle-resolving detector of that acceptance, or a slit). The beam
has a Gaussian spread σ_η in detuning centred on nq. Measured population:
P(t; σ_η) = ∫_channel p(δ) [I(nq + δ, t) − I(nq + δ, 0)] dδ / ∫_channel p(δ) dδ,
with the t = 0 (unmodulated) intensity subtracted ray by ray, which a measurement can do (it is
the pre-pump rocking curve). Reference: the flat-beam limit P_flat(t) = (1/q) ∫_channel [I(t) − I(0)] dδ.
Full TT `rocking_cube` on the channel grid once per case and geometry (spacing w_n/40, ≤ 1600 rays),
then every σ_η is a change of weights.

**Pre-registered rule (BACKLOG T2p, written before the run).** Four poster cases, Laue and Bragg,
α_ac = 20, step front, t̃ ≤ 1, σ_η = 0.3, 1, 3, 10 w_n. Metrics against P_flat: "period" = time of
the first maximum (Laue; parabolic interpolation) or the half-rise time (Bragg, monotone);
"amplitude" = P/P_flat at that time.
- **"The integrated observable works with an ordinary beam"** if at σ_η = 1, 3 and 10 w_n both
  |Δperiod| ≤ 1 % and |Δamplitude| ≤ 5 % in every case and geometry; otherwise the table of the
  smallest σ_η from which both hold (the integrated observable converges to P_flat from above in
  σ_η, the opposite direction to Entry 84), with w_n and q in µrad.
- Reported, not judged: P_flat's own shape against the Laue integrated-intensity law
  ∫_0^{2πα_n t̃} J_0(x) dx (normalised at t̃ = 1) and, for Bragg, against tanh²; and the signal size
  P_flat(1)/q relative to the centre intensity.
- Expectation written beforehand: at σ_η = 1 w_n the Gaussian still weights the satellite centre
  and P/P_flat will be above 1 by ~10–30 % with the centre law's period; from σ_η ≳ 3 w_n (the
  Gaussian flat over the satellite's ~2 w_n core) both criteria should hold; at 10 w_n the beam
  is 3–6 half-widths wide and is cut by the channel edge, which the normalisation handles. So the
  expected outcome is "condition table, σ_η ≳ 3 w_n", i.e. the integrated observable **does** lift
  the sub-µrad condition, with the detector acceptance (one channel, ~12 µrad) as the new condition.
- Convergence: δ-grid spacing halved at one case must change P/P_flat by < 0.2 %.

**Run.** `uv run python explore/nb/e86_population_divergence.py run fig` (43 s); data
`e86_population_divergence.json`, log `e86_population_divergence.log`, figure
`figures/paper/e86_population_divergence.png`. 359–1652 rays per case; δ-grid halving changes
P by 5e-5. Single carrier phase on the offset grid, as in Entries 1, 9, 84.

**Results under the pre-registered metrics** (Δperiod = first-maximum time in Laue, half-rise time
in Bragg, against the flat-beam channel integral P_flat; Δamp = P/P_flat at that time):

| geometry | (S, α, n) | w_n (half-widths) | σ_η = 0.3 w_n: Δperiod / Δamp | 1 w_n | 3 w_n | 10 w_n | smallest σ_η meeting both |
|---|---|---|---|---|---|---|---|
| Laue | (0.5, 2, 1) | 0.24 | +7 % / +940 % | −1.1 % / +790 % | −0.9 % / +360 % | −0.4 % / +65 % | none ≤ 10 w_n |
| Laue | (1.5, 2, 1) | 0.56 | +25 % / +580 % | +6 % / +380 % | −0.1 % / +113 % | 0.0 % / +12 % | none |
| Laue | (2, 4, 2) | 0.35 | +39 % / +360 % | +20 % / +250 % | +0.4 % / +65 % | +0.1 % / +6 % | none |
| Laue | (1.5, 4, 1) | 0.56 | +9 % / +290 % | −1.0 % / +160 % | −0.2 % / +28 % | 0.0 % / +2.4 % | 10 w_n |
| Bragg | (0.5, 2, 1) | 0.24 | +88 % / +440 % | +69 % / +420 % | +34 % / +280 % | +16 % / +57 % | none |
| Bragg | (1.5, 2, 1) | 0.56 | +65 % / +200 % | +43 % / +190 % | +10 % / +90 % | +1.4 % / +10 % | none |
| Bragg | (2, 4, 2) | 0.35 | +43 % / +150 % | +37 % / +130 % | +5.5 % / +56 % | +0.8 % / +5.8 % | none |
| Bragg | (1.5, 4, 1) | 0.56 | +82 % / +50 % | +69 % / +42 % | +9 % / +12 % | +0.9 % / +1.1 % | 10 w_n |

- **Verdict under the rule: condition table; the integrated observable does not pass at
  σ_η = 1, 3, 10 w_n.** Recorded as such, with the following caveat on the rule itself.
- **The amplitude metric was mis-specified in the pre-registration.** P is normalised per unit of
  beam inside the channel and P_flat per unit of channel width, so P/P_flat at a narrow beam is
  the ratio of the satellite's peak to the channel mean (the +940 % at σ_η = 0.3 w_n is the
  satellite peak being ten times the channel mean), i.e. it measures beam width over channel
  width, not the physics. A rule about the shape should compare curves normalised at the same
  time. This is a pre-registration flaw, stated, and the rule was not altered after the fact.
- **Post hoc, labelled:** the shape-normalised amplitude A = [P(t_ref)/P(1)]/[P_flat(t_ref)/P_flat(1)] − 1
  from the same traces, with the same 1 % / 5 % thresholds:

| geometry | (S, α, n) | 0.3 w_n: Δperiod / A | 1 w_n | 3 w_n | 10 w_n | smallest σ_η meeting both |
|---|---|---|---|---|---|---|
| Laue | (0.5, 2, 1) | +7 % / −36 % | −1.1 % / −14 % | −0.9 % / +1.7 % | −0.4 % / +1.7 % | 3 w_n |
| Laue | (1.5, 2, 1) | +25 % / +149 % | +6 % / +31 % | −0.1 % / +8 % | 0.0 % / +1.0 % | 10 w_n |
| Laue | (2, 4, 2) | +39 % / −25 % | +20 % / −3 % | +0.4 % / −0.2 % | +0.1 % / −0.1 % | 3 w_n |
| Laue | (1.5, 4, 1) | +9 % / +21 % | −1.0 % / +14 % | −0.2 % / +3.2 % | 0.0 % / +0.3 % | 3 w_n |
| Bragg | (0.5, 2, 1) | +88 % / −57 % | +69 % / −50 % | +34 % / −25 % | +16 % / −8 % | none |
| Bragg | (1.5, 2, 1) | +65 % / −52 % | +43 % / −42 % | +10 % / −15 % | +1.4 % / −2.3 % | none (10 w_n fails on period by 0.4 %) |
| Bragg | (2, 4, 2) | +43 % / −49 % | +37 % / −39 % | +5.5 % / −12 % | +0.8 % / −1.4 % | 10 w_n |
| Bragg | (1.5, 4, 1) | +82 % / −56 % | +69 % / −47 % | +9 % / −16 % | +0.9 % / −1.8 % | 10 w_n |

  In Laue the integrated population's *shape* is the flat-beam one from σ_η ≈ 3 w_n (three of
  four cases; 8 % at 3 w_n and 1 % at 10 w_n in the fourth): 3 w_n = 0.7–1.7 main Darwin
  half-widths = **1.3–3 µrad at Si 440 / 17.48 keV**, a ten-fold relaxation of the centre law's
  0.3 w_n and inside ordinary high-resolution beamline divergences, with the detector accepting
  one channel (q = 5–10 half-widths, 9–18 µrad there). In Bragg the satellite's rocking curve has
  1/δ² wings, the channel mean converges slowly, and even σ_η = 10 w_n (2.4–5.6 half-widths,
  i.e. a beam comparable to the channel) leaves 1–16 % on the half-rise time: the Bragg population
  is not a clean observable at any beam the channel can hold.
- P_flat's own shape is not the textbook law: against ∫_0^{2πα_n t̃} J_0 (Laue) the maximum
  deviation of the normalised curves is 13–25 %, against tanh² (Bragg) 25–38 %. The integrated
  population has a formation curve of its own (the channel integral of TT, which these runs
  tabulate) and no closed form here. Its signal is 5–29 % of the centre intensity per unit channel
  (P_flat(1) = 0.053–0.29).

**Answer to T2p.** Under the rule as written: no. In substance (post hoc, same traces): the
Laue integrated population is beam-tolerant from σ_η ≈ 3 w_n (1.3–3 µrad at Si 440 / 17.48 keV)
with a one-channel detector, which does lift the sub-µrad condition for the Laue train (P2); its
formation curve must then be taken from TT, not from sin² or ∫J_0. Bragg does not. A corrected
pre-registration (shape-normalised amplitude; Laue; σ_η = 3 w_n; the TT channel integral as the
law) is the item to run if P2's experimental E is to be raised; the data here are its first map,
not its test.

**Scope.** Step front, α_ac = 20, single phase, t̃ ≤ 1, Gaussian beam centred on the channel,
detector = one channel, t = 0 rocking curve subtracted ray by ray; lossless solvers; the ramp
regime and the ±η structure of the satellite not covered.

---

## 2026-10-23 — Entry 87: the Laue integrated population at σ_η = 3 w_n, corrected pre-registration (T2q)

**Why.** Entry 86's amplitude metric divided a beam-weighted mean by a channel mean, so it measured
the beam width, not the shape. BACKLOG T2q: the same question with the metrics written as formulas
first, normalisation denominators explicit. Script `explore/nb/e87_population_laue.py`.

**Definitions (fixed before the run).** Laue, α_ac = 20, q = 2πα_ac. Channel C_n = {η : |η − nq| ≤ q/2};
the detector accepts C_n and nothing else (q = 5–10 main Darwin half-widths for α = 4, 2, i.e.
9–18 µrad at Si 440 / 17.48 keV). Ray intensity change ΔI(η, t) = I(η, t) − I(η, 0⁺) (full TT
`rocking_cube`, the pre-pump rocking curve subtracted ray by ray). Beam p_σ(δ) ∝ exp(−δ²/2σ_η²),
δ = η − nq, σ_η = 3 w_n with w_n = 2πα|J_n(S)|.

    P_σ(t)  = ∫_{C_n} p_σ(δ) ΔI(nq + δ, t) dδ  /  ∫_{C_n} p_σ(δ) dδ          (beam-weighted population)
    P_0(t)  = (1/q) ∫_{C_n} ΔI(nq + δ, t) dδ                                    (flat-beam channel integral, the law)
    Ŝ_σ(t)  = P_σ(t) / P_σ(t_end),   Ŝ_0(t) = P_0(t) / P_0(t_end)                (shapes; denominator = value at the end of the trace, t_end = 1 t_tr for the step cases, 2 t_tr for the ramp)
    t*      = time of the first maximum of Ŝ (parabolic interpolation on the 1/60 grid)
    Δperiod = t*_σ / t*_0 − 1
    A       = Ŝ_σ(t*_0) / Ŝ_0(t*_0) − 1                                          (shape amplitude at the reference time)

**Pre-registered rule (BACKLOG T2q).** Four step cases (0.5, 2, 1), (1.5, 2, 1), (2, 4, 2),
(1.5, 4, 1) with t̃ ∈ [0, 1], plus one ramp case (1.5, 2, 1) with rise = 1 t_tr and t̃ ∈ [0, 2].
- **"The P2 integrated observable is measurable in the high-resolution beamline range
  (σ_η = 3 w_n = 1.3–3 µrad at Si 440 / 17.48 keV, one-channel detector 9–18 µrad)"** if
  |Δperiod| ≤ 1 % and |A| ≤ 5 % in all five runs; otherwise the condition table (σ_η = 3, 5, 10 w_n
  scanned for it, the smallest passing value per case).
- The law for the paper is the table P_0(t) itself (saved in the json), not sin² or ∫J_0.
- Expectation written beforehand, from Entry 86's post-hoc numbers: three step cases pass at
  3 w_n; (1.5, 2, 1) had A = +8 % there and may fail, passing near 5 w_n; the ramp is untested.
  So the likely outcome is "condition table with σ_η ≈ 3–5 w_n", not a clean pass.
- Convergence: δ spacing w_n/40 → w_n/80 at one case must change Ŝ_σ by < 0.2 %.

**Addendum to Entry 85 (physics cross-review, same day): the Bragg counterpart cannot use the
shear wave.** In symmetric Bragg 440 the reflection vector h is along the surface normal, and a
bulk wave travelling along the normal with transverse polarisation has u ⊥ h, so h·u = 0 and the
wave does not couple at all (the Laue poster geometry has h in-plane with u ∥ h, which is why it
couples). The Bragg window block was recomputed with a longitudinal [110] wave,
v_L = √((c11 + c12 + 2c44)/2ρ) = 9136 m/s (c11 = 165.7, c12 = 63.9, c44 = 79.6 GPa): f_D(Bragg)
= 0.423 GHz, x = 2.36 / 4.73 at 1 / 2 GHz, steady H = 0.059 / 0.123 (ratio 0.484), T_flat 2,
T_half 5, T_90 7 ns, 1/κ = 260 µm = 28 ns; the ringing-width column becomes 0.39 / 0.21 µrad
(1 / 2 GHz). Agreed with the reviewer; `e85_si440_case.py` and its log/json updated
(the Laue table is unchanged, since there the shear wave is the coupling one); the first
version's numbers are kept above, marked withdrawn. The elastic numbers for a longitudinal
[110] wave (ε = S d/Λ_ac with Λ_ac = v_L/f: 0.64× the shear strain at the same S and f; intensity
with ρv_L² in place of c44) are not re-tabulated here.

**Run.** `uv run python explore/nb/e87_population_laue.py run fig` (40 s); data
`e87_population_laue.json` (contains the law tables P_0(t) for all five runs), log
`e87_population_laue.log`, figure `figures/paper/e87_population_laue.png`. 359–1652 rays per run;
δ-spacing halving changes Ŝ_σ by 3e-6. Single carrier phase on the offset grid (the ramp's P_0
shows the carrier ripple; the metrics use the first maximum of the envelope-scale shape).

**Results** (Si 440 / 17.48 keV angles: 1 main Darwin half-width = 1.77 µrad):

| run (S, α, n, rise) | w_n (hw) | channel q | σ_η = 3 w_n: Δperiod / A | 5 w_n | 10 w_n | smallest passing σ_η |
|---|---|---|---|---|---|---|
| (0.5, 2, 1, step) | 0.24 | 10 hw = 18 µrad | −1.05 % / +2.7 % (period fails by 0.05 %) | −0.9 % / +2.9 % | −0.6 % / +1.9 % | 5 w_n = 2.1 µrad |
| (1.5, 2, 1, step) | 0.56 | 18 µrad | −0.07 % / +8.1 % | −0.02 % / +3.6 % | 0.0 % / +1.0 % | 5 w_n = 4.9 µrad |
| (2, 4, 2, step) | 0.35 | 5 hw = 9 µrad | +0.4 % / −0.2 % | +0.2 % / −0.2 % | +0.1 % / −0.1 % | 3 w_n = 1.9 µrad |
| (1.5, 4, 1, step) | 0.56 | 9 µrad | −0.2 % / +3.2 % | −0.1 % / +1.3 % | 0.0 % / +0.3 % | 3 w_n = 3.0 µrad |
| (1.5, 2, 1, ramp 1 t_tr) | 0.56 | 18 µrad | −0.6 % / +11.9 % | −0.2 % / +5.4 % (fails by 0.4 %) | 0.0 % / +1.5 % | 10 w_n = 9.9 µrad |

- **Verdict under the rule: condition table** (two step runs and the ramp fail at 3 w_n, as the
  pre-run expectation said for one of them). The table: every step run passes from σ_η = 5 w_n
  (2–5 µrad at Si 440 / 17.48 keV), the α = 4 runs already at 3 w_n (2–3 µrad); the 1 t_tr ramp
  needs 10 w_n (10 µrad, half the channel). The period is never the problem (≤ 1.05 % at 3 w_n,
  ≤ 0.9 % beyond); the amplitude A is, and it converges from above as 1/σ² (8.1 → 3.6 → 1.0 % for
  (1.5, 2, 1)), i.e. the residual weighting of the satellite core by the Gaussian.
- So the Laue integrated population is a beam-tolerant observable of the formation dynamics at
  **σ_η ≈ 5 w_n ≈ 2–5 µrad with a one-channel detector of 9–18 µrad** (step front), 10 µrad for a
  slow ramp: the range of an ordinary high-resolution beamline, 10–30× looser than the centre law
  (Entry 84) and 100× looser than the single-ray crystal factor (Entries 82–83). What is measured
  is not sin²(πα_n t̃): the flat-beam shape has its first maximum at α_n t̃ ≈ 0.4 with 1.6× the
  end value (step) and follows the TT channel integral P_0(t), tabulated in the json for the
  five runs; that table (or a TT run per case) is the law to compare an experiment with.
- The pre-registered threshold 3 w_n was too tight by a factor 1.7 for the α = 2 cases, where
  the channel (10 hw) holds the Gaussian comfortably but the satellite core is broader relative
  to q; the α = 4 cases pass. Not adjusted; recorded.

**Answer to T2q.** Not at 3 w_n for all five runs; yes at 5 w_n (step) and 10 w_n (ramp), which
is still within an ordinary high-resolution beam (2–10 µrad) and a 9–18 µrad detector channel.
For the paper: P2 as an integrated-population measurement in Laue is feasible with an ordinary
beam, with its formation curve supplied by TT; the E score of P2 can be raised to 4 on this
basis if Isaac accepts the 5 w_n condition (not raised by me, since the pre-registered 3 w_n
failed). Class: **confirmed, with a condition** (full TT beam average, five runs, metrics
defined in advance, no fit).

**Scope.** Laue, α_ac = 20, single phase, Gaussian beam centred on the channel, detector = one
channel, t = 0⁺ rocking curve subtracted ray by ray; lossless; n = 1, 2; a slow ramp of 1 t_tr
only. Bragg populations were shown unclean in Entry 86 and are not retested.

---

## 2026-10-23 — Entry 88: polarisation and the static Debye–Waller factor as rescalings (T2r)

**Why.** `notes/referee_pass.md` objection 2 leaves two model limits open: the π polarisation
component (DESIGN.md item 6) and the static Debye–Waller factor. BACKLOG T2r: do they reduce to a
redefinition of S and κ (α), or do they need numbers in the limitations section? Script
`explore/nb/e88_polarisation_dw.py`.

**Algebra written before the run.**
(a) Polarisation. The coupling carries C = 1 (σ) or |cos 2θ_B| (π); χ_0 does not. So for the π
component c → C c (both couplings; with χ_ih, c = iπα + κ_h/2 → C(iπα + κ_h/2)), κ unchanged:
α_π = C α_σ, g_π = g_σ/C, κ_h/κ unchanged, f_D,π = C f_D,σ so x_π = x_σ/C; the first-order
amplitude ν_π = C ν_σ, the formation rate α_n → C α_n. For a pure polarisation this is the
symmetric solver at another α: exact, as in Entry 81. An unpolarised beam adds the two
polarisation intensities, both driven by the same strain, so the oscillating parts add
deterministically: F_mix(ω) = F_σ(ω) + F_π(ω) (complex, each from CF-B with its own c), and the
mixed crystal factor H_mix = |F_σ + F_π| / (N_σ + N_π) with N the kinematic laws (N_π/N_σ ≈ C²
in the tail, since |R0| ∝ α there). For P2 the unpolarised centre intensity is
½[sin²(πα_n t̃) + sin²(πCα_n t̃)], two periods, not one law. C = 0.348 for Si 400 at 8 keV
(θ_B = 34.8°), 0.727 for Si 440 at 17.48 keV (21.7°), 0.866 for Si 400 at 17.48 keV (15.1°).
(b) Static Debye–Waller. χ_h = χ_h⁰ e^{−M} for both χ_rh and χ_ih (same factor), χ_0 unchanged:
α → α e^{−M}, κ unchanged, κ_h/κ ∝ e^{−M}. The X0h constants used throughout (Entry 20) are
already the room-temperature values with e^{−M} (Si 400: e^{−M} = 0.94 with B = 0.47 Å²), so the
model's Λ, g and χ_ih/χ_i0 already contain it; S = h·u is a phase and is untouched (the coherent
displacement enters through J_n(S), Entry 33). What is not in the constants is a temperature
change of M under drive: at the 1e-4 duty cycle of Entry 73 the average power is mW/cm² and the
change is negligible; stated, not computed.

**Pre-registered rule (BACKLOG T2r).** "Absorbed as a factor" if (i) the pure-π case is the σ
closed form with (g/C, x/C) to machine precision against the TT run at α = C (identity check at
one point), and (ii) the unpolarised H_mix from the closed-form sum agrees with the TT sum
(σ run + π run, same packet, Entry 75 machinery with c = iπC for π) to 0.3 % at every (f, T_w) for
Si 440 / 17.48 keV Bragg (longitudinal wave, x = 2.4–4.7, f = 1, 2 GHz, T_w = 7 and 300 ns) and
Si 400 / 8 keV (f = 20, 40, 100 GHz, T_w = 14 and 300 ns); otherwise "numbers in the limitations
section". Reported either way: H_π/H_σ, H_mix/H_σ, the window numbers T_flat/T_half/T_90 for σ,
π and unpolarised, and ν_π/ν_σ = C; for (b) the e^{−M} values and the statement above. No new
solver; the single-ray (σ-polarised) numbers of Entries 65–87 stand as the σ case.

**Run.** `uv run python explore/nb/e88_polarisation_dw.py` (3 s); data `e88_polarisation_dw.json`,
log `e88_polarisation_dw.log`. TT: Entry 81's `integrate` with c = iπ (σ) and c = iπC (π) on the
same packet and detunings, continuation past the packet as in Entry 75; closed form: CF-B with
the same two couplings; mixed = sum of the two complex Fourier components over the sum of the
two kinematic laws. (The window block's first version divided raw weights without the kinematic
normalisation and gave wrong steady ratios; fixed before the numbers below, rule unchanged.)

**(a) Polarisation.** C = 0.727 (Si 440 / 17.48 keV, θ_B = 21.7°), 0.349 (Si 400 / 8 keV, 34.8°).

| case | T_w | f | H_σ | H_π | H_π/H_σ | H_mix (unpolarised) | H_mix/H_σ | TT/CF (mix, π) |
|---|---|---|---|---|---|---|---|---|
| Si 440 / 17.48 keV, longitudinal | 7 ns | 1 GHz | 0.091 | 0.176 | 1.93 | 0.109 | 1.20 | 1.0000, 1.0000 |
| | | 2 GHz | 0.182 | 0.207 | 1.13 | 0.186 | 1.02 | |
| | 300 ns | 1 GHz | 0.059 | 0.115 | 1.93 | 0.078 | 1.32 | |
| | | 2 GHz | 0.123 | 0.229 | 1.87 | 0.159 | 1.30 | |
| Si 400 / 8 keV | 14 ns | 20 GHz | 0.911 | 0.950 | 1.04 | 0.911 | 1.000 | 0.9998, 0.9998 |
| | | 100 GHz | 0.949 | 0.950 | 1.00 | 0.949 | 1.000 | |
| | 300 ns | 20 GHz | 0.926 | 0.999 | 1.08 | 0.928 | 1.003 | |
| | | 100 GHz | 0.997 | 1.000 | 1.00 | 0.997 | 1.000 | |

Window numbers (CF-B, H(f_lo)/H(f_hi), targets 10 / 50 / 90 % of the way from 1 to steady):

| case | polarisation | steady ratio | T_flat / T_half / T_90 ns |
|---|---|---|---|
| Si 440 / 17.48 keV, 1 vs 2 GHz | σ | 0.484 | < 0.3 / 4.6 / 6.7 |
| | π | 0.501 | < 0.3 / 9.1 / 12.9 |
| | unpolarised | 0.491 | < 0.3 / 5.1 / 7.4 |
| Si 400 / 8 keV, 20 vs 100 GHz | σ | 0.929 | 5.0 / 12.9 / 24.6 |
| | π | 0.999 (flat: g_π = 0.13) | not defined |
| | unpolarised | 0.931 | 4.7 / 12.6 / 24.1 |

- **Rule (i), pure π:** identity to machine precision (TT/CF = 1.0000 at Si 440, 0.9998 at 8 keV,
  the latter the usual second-order residual at the packet amplitude used): the π component is
  the σ closed form with (g/C, x/C), ν_π = C ν_σ, formation rate C α_n. **Rule (ii), unpolarised:**
  closed-form sum against the TT sum, worst 0.021 % over all (f, T_w). **Verdict: absorbed as a
  factor** for a pure polarisation and as a two-term closed-form sum for a mixed beam.
- Size: at Si 400 / 8 keV the π component is 12 % of the weight (C² = 0.12) and H_mix is within
  0.3 % of H_σ at every window; at Si 440 / 17.48 keV it is 53 % and H_mix is 20–32 % above H_σ
  at 1–2 GHz, with window numbers 10 % longer than σ. So Si 400 / 8 keV can be measured with any
  polarisation; the Si 440 train needs the polarisation state known (a synchrotron or XFEL beam
  with a vertical diffraction plane is σ to better than 1 %, so this is a geometry choice, not a
  measurement). The formation law of an unpolarised beam is the two-period sum
  ½[sin²(πα_n t̃) + sin²(πCα_n t̃)], 0.3 % of a single law at 8 keV (C² weighting) and a visible
  second period at Si 440 (C = 0.73): same remedy.

**(b) Static Debye–Waller.** e^{−M} = 0.938 (Si 400), 0.880 (Si 440), 0.969 (Si 220) at 300 K
(B = 0.47 Å²). Both χ_rh and χ_ih carry it; χ_0 does not; so α → α e^{−M}, g → g e^{M},
χ_ih/χ_i0 → ×e^{−M}, and S = h·u is untouched. The X0h constants used since Entry 20 are the
room-temperature values with e^{−M} included, so every number in Entries 59–87 already has it;
nothing to re-tabulate. A change of M under drive heating at the 1e-4 duty cycle (mW/cm², Entry
73) is below any quoted precision; stated in the limitations, not computed. **Absorbed in the
constants.**

**Answer to T2r.** Both limits are absorbed: polarisation as a rescaling (α → Cα, κ fixed) with the
unpolarised case a two-term closed-form sum verified against TT to 0.02 %, Debye–Waller already
inside the X0h constants. The limitations section keeps one sentence each (polarisation state must
be known where C² is not small; DW is the room-temperature value). `notes/referee_pass.md`
objection 2 is now answered on all its items except the lateral acoustic uniformity (DESIGN.md,
not modelled). Class: **confirmed** (algebra + TT identity + 0.02 % mixed check).

**Scope.** Symmetric Bragg, first order, centrosymmetric; the mixed-beam sum assumes the two
polarisations see the same strain and detuning (true for a single ray; a divergent beam is
Entry 82's question, independent of polarisation). The P2 (Laue) two-period statement is algebra
from the rescaling, not run.
