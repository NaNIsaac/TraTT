# Paper draft skeleton (BACKLOG T2h, draft tier; nothing here is new computation)

Written 2026-10-23 from `notes/paper_candidates.md` and lab notebook Entries 1–81. Every number
cites the entry and the script that produced it; a claim's status is the notebook's (confirmed =
full TT in this model, converged, no fit; prior art = in literature read in full; hypothesis /
predicted = stated or computed, not tested against experiment). Posting, journal choice and
authorship are Isaac's. Prior-art items marked **unread** are on BACKLOG T2d.

## 0. Title candidates

1. Time-resolved Bragg reflectivity of an acoustic strain field: a windowed transfer function and
   the dynamical formation of satellites
2. What a time-resolved Bragg reflection measures of a GHz strain field: crystal factor, window
   condition and satellite formation
3. The crystal factor of time-domain X-ray phonon spectroscopy, from picosecond packets to
   nanosecond trains

(1) is the P1 + P2 paper; (3) is P1 alone. The choice follows Isaac's decision in
`paper_candidates.md` §4.

## 1. Abstract (150 words, P1 + P2 version)

An acoustic wave in a Bragg-diffracting crystal produces satellites that are dynamical reflections. Solving the two-beam Takagi–Taupin equations for a strain field entering a
semi-infinite crystal, we show to first order in the strain that the time-Fourier weight of the
reflectivity at the phonon frequency is the strain spectrum times a closed-form crystal factor
H(ω; T_w), of the measurement window T_w, the absorption parameter χ_i0/|χ_rh| and the anomalous
ratio χ_ih/χ_i0, with nothing fitted. H is flat for windows shorter than the strain's transit of the
extinction and absorption depths and steady beyond; for Si 400 at 17.5 keV the crossover is 20–50 ns.
Anomalous absorption makes the two detunings of one phonon frequency unequal, a ratio free of the
strain spectrum. Beyond first order the satellites form with extinction length Λ_P/|J_n(S)|, and the
transient reads a one-transit running mean of J_n(S A). We give the window and amplitude conditions for each regime.

(Count: 152 words. The P1-only abstract drops the last two sentences and adds the ps-packet delay
scan.)

## 2. Sections

1. Introduction: time-resolved X-ray diffraction of acoustic strain (ps ultrasonics, bulk GHz
   waves); what is read from the time-domain signal; prior art list (§6); the three questions
   (crystal factor, window, amplitude).
2. Model: two-beam Takagi–Taupin, Bragg Riccati equation in retarded time (Entry 66 trick),
   normal and anomalous absorption (η → η + iκ; c → iπα + κ_h/2, Entry 77), asymmetric geometry
   as a rescaling (Entry 81), validity (quasi-static, plane wave, σ, centrosymmetric). DESIGN.md
   items 1–7 restated as the assumptions.
3. First order: the windowed closed form. 3.1 derivation (Entry 59 linearisation, Entry 75
   integration by parts, CF-B); 3.2 where the plane-wave form fails (Entry 64 vs Sondhauss–Wark);
   3.3 the window condition and the three design numbers T_flat, T_half, T_90 (Entry 75); 3.4 the
   crystal parameters g and χ_ih/χ_i0 (Entries 65, 77); 3.5 the ±η side asymmetry (Entries 77,
   78); 3.6 the linearity bound ν ≲ 0.1 on its domain (Entry 68); 3.7 **the beam condition**: the crystal
   factor is a single-ray property; a beam spread σ_η ≫ κ + π/x dephases the ringing after 1/σ_η and
   leaves the kinematic law times a flat κ√(π/2)/σ_η at every window (Entry 82, TT beam average vs
   closed form 0.02 %), so H(f; T_w), χ_ih and the ±η asymmetry need ≲ 0.03–0.4 µrad and ≲ 2–4 meV.
4. Beyond first order: formation law I_n = sin²/tanh²(πα|J_n|t̃) (Entries 1, 9, 19, 24), the
   running-mean reading Θ_n and its degeneracy class with the R ≤ 0.1 boundary (Entries 28–34,
   72, 76, 79), the ramp regime and overshoot map (Entry 67). [P2; a section of P1 if P1 alone.]
5. Experimental design (proposal, not data). 5.1 picosecond packet with a delay scan ≥ T_90
   (P1 alone; Reis-type geometry, Entry 75: H(20)/H(100) = 0.558 at 52 ns for a 10 ps packet);
   5.2 GHz bulk shear train in Si 440 Laue / Si 400 Bragg read over tens of ns (P1 + P2, Entries
   20, 73, 80: S ≤ 6 clean on the elastic side, drive power 6–900 W/cm² peak unverified);
   5.3 **asymmetric cut to shorten the delay scan**: b = γ_0/|γ_h| = 4 brings T_90 from 52 to
   33 ns, b = 10 to 16.5 ns, at the cost of the γ_0/|γ_h| intensity factor and a 2.7° exit
   (Entry 81); 5.4 reading χ_ih/χ_i0 from a ±ω satellite pair, sensitivity 0.02 (17.5 keV) to
   0.07 (8 keV) per unit, |R0| ratio 1.0000 in the phonon band (Entry 78); 5.5 what to measure
   and to what precision (window known to 10 % on the steep part, Entry 75); 5.6 **beam resolution**:
   without meV / sub-µrad resolution the measurement is the flat-factor regime of Entry 82 and the
   estimator needs no crystal factor; with it, 5.1–5.4 apply. State which regime the proposal is in.
6. Discussion: relation to Sondhauss–Wark 2003 (Entries 64, 69), to Chalise et al. (Entries 57,
   66, 71: their 140 ps window sees no crystal factor; the second-order term K = 0.65–0.88
   remains), what is not new (Bessel ladder, Pendellösung), what is (window condition, two-
   parameter H verified inside the Darwin width, degeneracy class, side asymmetry), and the
   prior-art caveat (§6, unread items).
7. Limitations (§5 below).
Appendix A: the retarded-time Riccati and the analytic continuation past the packet (Entry 75).
Appendix B: pre-registered rules and their outcomes, one line each (Entries 68, 70–72, 75–79, 81).

## 3. Figures 1–4, data sources

| figure | panels | source (entry; script; reproduce) | status |
|---|---|---|---|
| Fig. 1 estimator and one-curve collapse | (a) R(t) at one detuning for a packet with the Fourier component at ω = \|η\| marked; (b) W/W_kin against ω/ω_D for three packet spectra at κ/2πα = 0.12 | Entry 59; `explore/nb/paper_figs.py fig1` → `figures/paper/fig1_estimator.png` (a: e75 packet run, 0.2 s; b: e58 json + e59 closed form); captions in `figures/paper/CAPTIONS.md` | confirmed |
| Fig. 2 where the plane-wave first order fails | \|F_pred\|/\|F_TT\| against ω/ω_D for S&W eq. 29, S&W normalised, Bloch-wave (Riccati) form; near-peak error | Entries 64, 69; `explore/nb/paper_figs.py fig2` → `figures/paper/fig2_planewave_vs_bloch.png` (from `e64_sw.json`) | confirmed in our variables; S&W PDF not seen (Entry 69, transcription of Entry 64) |
| Fig. 3 the window and the crystal parameters | (a) H(f; T_w) against T_w with TT circles and T_flat/T_half/T_90; (b) H(20)/H(100) against T_w; (c) H(f) with and without χ_ih at both energies; (d) ±η asymmetry against f at four windows | (a, b) Entry 75; `explore/nb/e75_window_cf.py run posthoc fig` → `figures/paper/e75_window_H.png`. (c) Entry 77; `e77_chi_ih.py run ttp ttm conv60 fig` → `figures/paper/e77_chi_ih.png`. (d) Entry 78; `e78_side_asymmetry.py run fig` → `figures/paper/e78_side_asymmetry.png` | confirmed (a–c), model-confirmed prediction (d) |
| Fig. 4 linearity bound and the dynamical satellite | (a) 1 − K against ν for the 24 series of Entry 68 with the tanh² and sin² laws; (b) the formation-law collapse I_n(t̃) = sin²(πα\|J_n\|t̃) | (a) Entry 68; `e68_nu_collapse.py run stats fit2 fig` → `figures/paper/e68_nu_collapse.png`. (b) Entries 1, 9; `explore/nb/paper_figs.py fig4b` → `figures/paper/fig4b_formation_collapse.png` (from `explore/nb/poster_f1_points.npz`; poster frozen) | confirmed |

Supplementary candidates: `figures/paper/e72_history_degeneracy.png`, `e76_degeneracy_S.png`,
`e79_degeneracy_R.png` (degeneracy class and its R boundary); `figures/si400_harmonic.png`
(second-order term for a ps packet, Entry 66); `figures/poster/f7_overshoot.png` (Entry 67, frozen).

## 4. Claims and their status

| # | claim | status | evidence |
|---|---|---|---|
| 1 | First-order linearisation of TT for a phonon | prior art | Sondhauss & Wark 2003 eq. 25–29 (Entries 64, 69) |
| 2 | Closed-form H with the Bloch-wave (Riccati) depth dependence, valid inside the Darwin width where the plane-wave form is off by 3–8× | confirmed in the model; novelty bounded by **unread** Chukhovskii & Förster 1995, Wark & Lee 1999 | Entries 59, 64, 65 |
| 3 | H needs a long window: the strain must travel ≥ max(1/κ, x/π) Λ during it; closed form H(f; T_w) exact to first order at every window 0.14–300 ns | confirmed (TT 0.3 %); packet-free short form rejected below ~100 τ | Entries 71, 75 |
| 4 | Design numbers Si 400 / 17.48 keV: T_flat 19.6, T_half 32.8, T_90 52.2 ns; ripple with period 2x Λ/v | confirmed (TT to 0.002) | Entry 75 |
| 5 | Chalise et al. 140 ps window sees no crystal factor (H flat to 4 %); their low-frequency deficit is only partly the second-order term K = 0.65–0.88 | confirmed in the model; their normalisation unknown; any real beam also flattens the crystal factor (Entry 82) | Entries 66, 71, 82 |
| 6 | χ_ih is a second tabulated crystal parameter: 0.9 % steady, 1.4 % at 52 ns, 3–4 % at 8 keV | confirmed; pre-registered rule fired | Entry 77 |
| 7 | ±η side asymmetry, strain-spectrum-free, reads χ_ih/χ_i0 (1.9–6.2 %, sensitivity 0.02–0.07) | predicted and confirmed in the model (TT 0.06 %); unobserved; a single-ray observable: odd in δ, it averages to zero over a symmetric beam (Entry 82) | Entries 77, 78, 82 |
| 8 | Asymmetric Bragg = exact rescaling; one geometry factor g(√b + 1/√b)/2; T_90 shortens with Λ_asym | confirmed (TT 1e-16) | Entry 81 |
| 8b | Beam divergence / bandwidth: for σ_η ≫ κ + π/x the crystal factor is flat, κ√(π/2)/σ_η, at every window (x ≳ 10); H(f; T_w) and the ±η asymmetry are single-ray (meV, sub-µrad) observables; at x = 1–5 the tolerance is 0.3–0.6 µrad / 4–7 meV (8 keV) and the beam-averaged factor is non-monotonic in σ | confirmed (TT beam average vs closed form 0.02–0.24 %) | Entries 82, 83 |
| 9 | ν ≲ 0.1 ⇒ weight linear to 2 % (1 − K ≤ 1.8 ν²) | confirmed as a bound on the scanned domain (Bragg, x ∈ {2, 5}, one packet family, weight at the peak) | Entry 68 |
| 10 | "One dimensionless criterion joins the regimes" | falsified (spread ×10–14 at fixed ν) | Entry 68 |
| 11 | Formation law I_n = sin²/tanh²(πα\|J_n\|t̃), Bessel ladder, period-ratio calibration | confirmed; prior art in kind DeCamp 2001 (n = 0), Polikarpov 1987/88 (static); as a centre observable it needs the beam condition of 12b | Entries 1, 9, 13, 19, 24; tests |
| 12 | Order n reads the one-transit running mean of J_n(S A); zero-mean perturbations degenerate (≤ 3e-6 for S ≤ 6); boundary R = 0.2 | confirmed with a boundary | Entries 72, 76, 79 |
| 12b | Beam condition for the centre formation law: σ_η ≲ 0.3 w_n (0.1–0.3 µrad, Si 440 / 17.48 keV Laue; 0.6–1.4 µrad, Si 400 / 8 keV Bragg); at σ_η = w_n the first maximum loses 13–40 % | confirmed (full TT beam average, 8 cases) | Entry 84 |
| 8c | Polarisation: π = σ form with α → Cα (κ fixed), unpolarised = two-term closed-form sum (TT 0.02 %); negligible at Si 400 / 8 keV (C² = 0.12), 20–32 % at Si 440 / 17.48 keV; static Debye–Waller already in the X0h constants | confirmed | Entry 88 |
| 12c | Integrated Laue satellite population: flat-beam shape from σ_η ≈ 5 w_n (2–5 µrad, step) and 10 w_n (10 µrad, 1 t_tr ramp) with a one-channel detector (9–18 µrad); formation curve from TT (tabulated; 13–25 % off ∫J_0); Bragg population not clean | pre-registered threshold falsified: 3 w_n failed in 3 of 5 runs (Entry 87; Entry 86 failed under its own rule). 5 w_n (step) and 10 w_n (ramp) are post-hoc conditions from the same cases; confirmation needs a pre-registered run on new cases (T2s) | Entries 86, 87 |
| 13 | Signal maximum = Fisher-information maximum for a ps packet at one angle | negative result | Entry 70 |
| 14 | Bessel-zero resurrection at third order | confirmed in the model, unobservable (4e-5); **unread** Entin, Iolin–Entin | Entries 45, 51, 52, 54 (not in this paper; P4) |
| 15 | Realisability: S ≤ 6 at 1–2 GHz clean on the elastic side; drive power 6–900 W/cm² peak | constants only; harmonic column withdrawn (shear symmetry); transducer number unverified | Entries 73, 80 |

## 5. Limitations (to be stated in the paper)

- Two-beam, plane-wave, quasi-static (β_qs ≤ 1e-2), centrosymmetric crystal (χ_h = χ_hbar); σ
  polarisation throughout, with π as the rescaling α → Cα and an unpolarised beam as the two-term sum
  (Entry 88: state the polarisation where C² is not small); Debye–Waller at room temperature inside
  the X0h constants; beam divergence and bandwidth as a single-ray average (Entries 82–87).
- Bragg: symmetric or asymmetric by rescaling, normal and anomalous absorption; Laue: symmetric,
  lossless only (χ_h'' excluded in Laue; DESIGN.md item 4).
- First order in S for the closed form (ν ≪ 0.1); the ±η ratio is second-order sensitive and was
  verified at peak strain 1e-6 (Entry 78).
- The strain enters only as a phase h·u (no elasto-optic change of χ_h, no Debye–Waller
  reduction); packet enters from the surface at t = 0; windows > 14 ns put the front 0.1–2.5 mm
  deep, so the crystal must be that thick.
- Inside the Darwin width (x ≲ a few) the |R0| asymmetry (10 %) must be carried in the ±η ratio;
  all phonon-band numbers are at x ≥ 27.
- All single-ray results (H(f; T_w), χ_ih, ±η asymmetry, design numbers) hold for a beam spread in
  detuning below κ + π/x; for an ordinary beam the crystal factor is flat (Entry 82).
- The centre formation law (Fig. 4b) holds for σ_η ≲ 0.3 w_n (Entry 84); an ordinary beam measures
  the integrated satellite population instead, whose law is not given here.
- The Chalise comparison is "same sign, same order": their energy, window and normalisation are
  not published in the form needed.
- Experimental: no long-window satellite measurement exists; transducer peak power unverified;
  the S = 3–6 regime is limited by drive power only (Entry 80), not by the wave.
- Numerical: every entry carries a convergence check (step, phases, S); the 0.3 % residual of the
  first-order comparisons is the second-order term at S/100.

## 6. Prior art (bibliography to be checked before citing; "unread" = not read in full here)

- Sondhauss & Wark, Acta Cryst. A 59, 7 (2003): first-order TT for a phonon; plane-wave form fails
  within ~1 Darwin width (Entries 64, 69). Read from a transcription; PDF not seen (human-todo).
- Chukhovskii & Förster 1995; Wark & Lee 1999: time-dependent TT, possible earlier Bloch-wave
  first-order form. **unread** (DotFiles#55).
- Chalise et al. (Entry 57): ps-ultrasonics strain spectrum vs reflectivity spectrum, Fig. 2.
- Shayduk 2013: names the heterodyne oscillation (Entry 60).
- DeCamp et al. 2001 (n = 0 pulse); Polikarpov 1987/88 (static extinction length); Bucksbaum &
  Merlin 1999 "phonon Bragg switch" **unread**; Zolotoyabko, Iolin on acoustically modulated
  crystals **unread** (exact papers to identify); Grigoriev et al. APL 2006 piezoelectric X-ray
  switch **unread**.
- Reis et al. PRL 86, 3072 (2001) **unread**; Trigo et al. Nat. Phys. 9, 790 (2013) **unread**;
  Schick et al. udkm1Dsim CPC 185, 651 (2014) **unread**: the P1 comparisons a referee will name.
- Entin JETP Lett. 26, 392 (1977), JETP 77, 214 (1979); Iolin & Entin JETP 85, 1692 (1983):
  X-ray–acoustic resonance, overlap regime **unread** (P4 only).
- Zarembo & Krasil'nikov, Sov. Phys. Usp. 13, 778 (1971): no second harmonic of shear waves
  **unread** (Entry 80).

## 7. What is still missing before a first full draft

- Isaac's decision P1 + P2 vs P1 alone (title, §4 and §5.2 depend on it).
- (done, T2i) Fig. 1, Fig. 2, Fig. 4b are files in `figures/paper/` with captions in `CAPTIONS.md`.
- The three unread 1990s/2001 papers for the novelty sentences of §3 and §6.
- A transducer number for §5.2 (DotFiles#57).
