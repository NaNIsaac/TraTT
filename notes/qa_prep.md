# Poster session prep — P2.17, KPS Fall 2026, Jeju, 22 October

Every answer below carries the entry in `lab_notebook.md` that backs it. If a question has no
entry, the honest answer is "not checked" — say that, and write the question down.

Rule for the session: **concede the known part first, then say what is added.** Every strong
question here has a real prior-art or validity component. Conceding it costs nothing and makes
the rest credible.

---

## 1. The opening

**10 seconds.** "Acoustic satellites are known in steady state. We resolve how they *form* while
the acoustic wavefront is still crossing the crystal — and each order turns out to have its own
clock."

**30 seconds (Korean, for a first approach).** "음파가 결정을 지나는 동안 위성 피크가 어떻게
생기는지를 계산했습니다. 핵심은 파면 뒤쪽이 '채워진 깊이짜리 결정'처럼 행동한다는 것입니다.
차수 n의 결합이 πα|J_n(S)|이라서, 차수마다 형성 속도가 다르고 그 순서를 음향 진폭 S로
바꿀 수 있습니다. 맞춤 변수는 없습니다."

**60 seconds, pointing at panel 3.** "Behind the front, satellite *n* is a crystal of the filled
depth v_s t whose coupling is πα|J_n(S)|. So the Laue centre intensity is sin²(πα_n t) and Bragg
is tanh²(πα_n t), with α_n = α|J_n(S)|. Four different (S, α, n) collapse onto one curve, no fit
parameter. The transit time doesn't set the formation time — it only ends the transient. The
clock is t_tr/α_n, the satellite's own extinction length." {E1}{E9}{E33}

---

## 2. Prior art — the questions that decide whether you are taken seriously

### Q. "Pendellösung in the time domain is not new. DeCamp did this in 2001."

**Concede, precisely.** Yes — DeCamp, Reis & Bucksbaum, *Nature* **413**, 825 (2001), read in
full: Ge (001), 280 µm, asymmetric 202 Laue at 10 keV, a **single short strain pulse** acting as
a moving interface between two perfect crystals, switching the **main** reflection. Their own
statement is that the pulse must travel about a quarter of a Pendellösung depth. {E21}

**Then.** "No satellites, no periodic drive, no Bessel functions in that work. Our n = 0 law is
the CW-grating version of their picture and we cite it on the poster. The satellite results sit
outside it: every order n, both geometries, one law." {E21}{E7}

### Q. "The satellite extinction length ~ Λ/|J_n| is in the ultrasonic-superlattice literature."

**Concede fully.** Correct — Polikarpov, *phys. stat. sol. (b)* **143** (1987) and **149** (1988),
read in full; Iolin & Entin; Zolotoyabko and co-workers. That is the steady-state ingredient and
we attribute it. {E29}{E50}{E7}

**Then.** "What we add is the transient: the time-domain collapse onto the static thickness
curves, formation ordered by the Bessel ladder, the accumulated-coupling law for a finite ramp,
and the path memory." {E21}

### Q. "Zolotoyabko already did stroboscopic X-ray diffraction from SAWs."

Standing SAW at 0.58 GHz in LiNbO₃, satellites modulated at 2f — a phase-locked **periodic steady
state**, not a formation transient. Different object. {E21}

### Q. "Why not do this with a surface acoustic wave? That is how everyone does it."

Expect this one. Say the geometry first, then the number. {E63}

**Geometry.** "A SAW's wavevector lies *in* the surface and its amplitude decays into the depth.
Ours is a bulk wave travelling *into* the depth, uniform behind the front. So what grows with time
is different: for us the modulated depth, for a SAW the modulated lateral extent. Our solver is
one-dimensional in depth and cannot represent an in-plane wavevector at all — which is why that
literature talks about sagittal diffraction."

**The number that settles it.** "A Rayleigh wave's skin is 0.30–0.41 Λ_ac, so the fraction of the
extinction depth it modulates is 0.30–0.41 R — independent of frequency and material. Resolving
the channels needs R ≤ 0.34–0.57, so the modulated fraction is capped at 10–17 %. A SAW cannot
fill the extinction depth and keep its satellites resolved: resolving wants a short acoustic
wavelength, filling wants a long one."

**The physical point, if they push.** "The skin depth is set by the frequency and does not grow.
Once the surface wavefront has passed a point, the profile under it is static. **A SAW has no
depth-filling transient** — that is the whole mechanism here. A SAW does have a transient, the
lateral filling of the X-ray footprint, ~20 ns for a 100 µm footprint. Same order as our t_tr,
but it is an inhomogeneously illuminated surface, not a depth-graded crystal."

**What does carry over.** "The accumulated-coupling law is a statement about the depth *profile*,
not about the front. I checked it on a static exponential skin: it holds to ≤ 1.1×10⁻² absolute
for skins thinner than 0.1 L, which is the whole SAW regime. And a skin of depth d reads like a
filled depth of 1.0–3.1 d depending on S — above S ≈ 1.8 the decaying tail sweeps through the
maximum of J₁, so depth where the amplitude has *fallen* couples more strongly than the surface."

**Close by naming your own geometry.** "So the experiment this predicts is a bulk wave from a
transducer, launched along the depth, and shear-polarised — in symmetric Laue a longitudinal wave
along the normal gives h·u = 0." {E6}{E20}

Never say a SAW experiment is uninteresting. It measures a different quantity.

### Q. "What about Sondhauss & Wark's time-dependent dynamical diffraction theory?"

**This one only touches the Outlook line, so route it there.** "Acta Cryst. A59 (2003) — closed
access, requested through the library, not read. Its abstract says it contains perturbative
analytical solutions, so I do not claim the first-order solution step. What has no visible
precedent is using it as a transfer function: the field's own 2023 review doesn't cite the paper
and contains no occurrence of Takagi, Darwin or extinction, and the 2026 paper that would need
the correction doesn't have it." {E60}{E61}

### Q. "So what is actually yours?"

Say exactly these five, nothing more: {E21}

1. transient formation of satellites behind a travelling front;
2. the time-domain collapse onto the static thickness curves, no fit parameter;
3. formation order set by the Bessel ladder and reorderable with S;
4. the accumulated-coupling law Θ_n for a finite ramp;
5. the period-ratio calibration of S.

Add: "The search was web-level, not a database search. I would like to know if you have seen any
of these." — turn it into an ask.

---

## 3. The model — the sharpest technical question

### Q. "Where is the time derivative of the X-ray field? This is not time-dependent diffraction."

**Correct, and it is deliberate.** The deformation is frozen at each instant and the two-beam
Takagi–Taupin system is solved in depth. The number: in the poster's window the X-ray crosses the
crystal in 0.99 ps, which is **0.2 % of a 500 ps acoustic period**; β = 2πfL/c ≤ 1.2×10⁻² over
Si 440, 17–30 keV, L = 163–296 µm, 0.6–2 GHz. The lattice is static while the field is inside.
{E62}

**What it does not cover, said before they say it.** The crystal's own field settling time — the
Shvyd'ko–Lindberg regime — is a separate question, and this model does not address it. Any
dependence on v_s beyond t_tr is outside these maps by construction. {E21}{E62}

### Q. "Are these Floquet sidebands / Stokes and anti-Stokes?"

"The channels sit at η = nq in the rocking variable. There is no photon-energy axis in this
calculation, so I don't label them Stokes / anti-Stokes — that needs the frequency-resolved
treatment." Do not improvise past this. {E0}

### Q. "Absorption?"

Simulated, not estimated: `core.solve_laue` takes a complex α. The oscillation contrast goes as
1/cosh(μLκ|J_n|t). Si 220 at 8 keV with L = 100 µm gives μL ≈ 1.4 and contrast ≈ 0.78 at t = t_tr.
The recommended window is better — Si 440 at 17–30 keV has α″/α′ < 0.4 % and contrast ≈ 1.00.
Design rule: keep μLκ|J_n| ≲ 1. {E11}{E20}

### Q. "How well does the law hold, and where does it break?"

Panel 6 is that map, and the limit is **R = Λ_ac/Λ_P**, not α: {E26}{E53}{E54}{E25}

- accuracy ≈ 0.9 R, independent of how many turns the clock makes — f tunes accuracy, L and S
  tune the clock;
- with the exact acoustic-phase average, channels stay independent to better than 5 % up to
  R = 0.34–0.57 depending on S;
- the law fails when |J_n| α_ac ≲ 0.5 — at a Bessel zero the slow clock stops and only an
  acoustic-period ripple ~ (α/α_ac)² is left;
- near R ~ 1 a Bessel-dark order reappears through a third-order coupling πα(R/2)²C_n(S), with a
  renormalised zero. {E45}{E52}

### Q. "Is it a special-point artefact — integer α_ac, nice α?"

Checked on random conditions, not a grid. {E19}

### Q. "How is the solver verified?"

54 tests. Analytic limits to 10⁻⁹ (whole rocking curve against the closed form, Bragg tanh²(πα),
RK4 fourth-order convergence measured at 3.98–4.00). An independent exact Floquet solver agrees to
3×10⁻⁷; unitarity holds to 3×10⁻⁷. {E0}{E42}

---

## 4. Numbers a visitor may challenge

| They ask | Answer | Entry |
|---|---|---|
| Why 0.95 t_tr/α_n? | 99 % formation of tanh², 0.953/α_n, capped at t_tr — it cannot exceed the transit | E9 |
| Why does S reorder the channels? | t\* / t_tr = 1/(2α\|J_n(S)\|); at S = 3 that is n = 2, then 1, then 3 | E1 |
| The dark point — is the satellite gone? | Centre drops to 0.1–0.2 %, but ~52 % of the integrated satellite is still there, in two side lobes; lobes collapse as (η−nq)t = G(α_n t) | E23 |
| Centre vs integrated intensity? | Same period, different phase and damping: the integral leads by ~1/8 period and damps as J₀ | E23 |
| Is the transient-only satellite just a slow amplitude sweep? | See §5 — concede most of it | E34 |
| 1–3 % on the phase walls? | Formula-level check against the plane-wave stacking-fault result with δ = n·dφ (Authier & Simon 1968; Manukyan 2014) | E56 |

---

## 5. The one claim to state carefully

**"Transient-only satellite: an order forbidden at the final amplitude is the brightest during
formation."** The honest framing is already written down — use it verbatim. {E34}

- For a **slow** ramp (T ≫ t_tr) this is the steady Bessel ladder sampled on the way to a zero.
  Any amplitude sweep shows n = 1 at intermediate S. Say this first.
- What the transient adds: (i) the t_tr memory — the satellite reports the drive t_tr/2 late and
  smoothed; (ii) the fast-ramp regime T ≲ t_tr, where it is a flash set by the ramp profile
  *inside* the crystal; (iii) a step onset gives **exactly zero** at all times, so the peak is a
  signature of the amplitude gradient inside the crystal, not of the final state.
- Robustness: > 30:1 transient-to-final contrast survives a 5 % error in S and 10 % acoustic loss
  across the crystal.

---

## 6. Experiment — for the person who might actually do it

Verified X-ray constants from X0h; c44 and the transducer are still estimates. {E20}

| | Si 440, 17.48 keV | Si 440, 30 keV |
|---|---|---|
| Λ_P | 54.3 µm | 98.8 µm |
| L = 3Λ_P | 163 µm | 296 µm |
| t_tr | 27.9 ns | 50.7 ns |
| τ₁ | 16.7 ns | 30.3 ns |
| f_min (r = 0.1) | 1.1 GHz | 0.6 GHz |
| strain | 2.7×10⁻⁵ | 1.5×10⁻⁵ |
| transmission | 0.78 | 0.92 |

Points to make: tens of ns is comfortable against ~100 ps storage-ring bunches; higher energy and
higher order improve **every** axis at once (lower frequency, lower strain, less absorption).

**"How do I know S in the experiment?"** From the ratio of the two temporal periods,
τ₂/τ₁ = |J₁/J₂| — independent of gain and of thickness, ~1 % recovered in a simulated measurement
with 3 % noise. The design rule is α|J₂(S)| ≥ 0.5; the S = 1 case failed at +9.4 % precisely
because α₂ = 0.35 and the n = 2 trace never completed half an oscillation. {E13}

---

## 7. Negative results — volunteer these, they buy credibility

- The satellite **peak shift** is kinematic interference, not a dynamical shift. {E3}
- The transient **±n asymmetry** is a first-acoustic-wavelength effect ∝ 1/α_ac, not dynamical.
  {E8}
- Same strain **up vs down** a pulse: mostly a negative result for the whole rocking curve. {E36}
- **Chirped** acoustic waves give a Landau–Zener plateau, no focusing. {E38}
- The session-1 conclusion "no additional timescale" was **wrong** — t_tr/α_n is one. Saying this
  out loud is a strength. {E12}

---

## 8. Do not say

| Never | Say instead |
|---|---|
| "time-dependent Takagi–Taupin" | "Takagi–Taupin with a travelling deformation, solved at each instant" |
| "Stokes / anti-Stokes", "Floquet sidebands" (energy sense) | "channels at η = nq"; no photon-energy axis here |
| "we derived the transfer function" (as new) | "a closed-form transfer function — prior art not yet settled, Sondhauss & Wark requested" |
| "temporal Pendellösung is new" | "the organising principle, from DeCamp 2001, extended to every order" |
| "satellite extinction length ~ Λ/\|J_n\| is ours" | "that is the steady-state ingredient, from Polikarpov and the superlattice literature" |
| any number not on this sheet | "I don't have that checked — can I take your email?" |

---

## 9. Pocket card

```
clock        t_tr / alpha_n,   alpha_n = alpha |J_n(S)|
Laue         sin^2(pi alpha_n t)        Bragg  tanh^2(pi alpha_n t)
order        t*/t_tr = 1 / (2 alpha |J_n(S)|)      S = 3 -> n = 2, 1, 3
dark point   centre 0.1-0.2 %, integral still 52 %
validity     R = Lambda_ac / Lambda_P <= 0.34-0.57;  accuracy ~ 0.9 R
verified     1e-9 analytic, 3e-7 Floquet, 54 tests
experiment   Si 440, 17-30 keV, 1-2 GHz, tau_1 = 17-30 ns, strain ~3e-5
calibrate S  tau_2 / tau_1 = |J_1 / J_2|,  needs alpha |J_2| >= 0.5
concede      DeCamp 2001 (main line, pulse);  Polikarpov / Entin / Zolotoyabko (steady state)
```
