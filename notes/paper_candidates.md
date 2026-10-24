# Paper candidates — what in TraTT can carry a journal paper

Written 2026-10-22 (night), for Isaac's directive of the same day (`docs/isaac-directive-20261022.md`
in the NaNIsaac workspace), BACKLOG T1d. Every statement points to a notebook entry. Status words:
**confirmed** = full TT in this model, converged, no fit; **prior art** = in the literature we have
read in full; **hypothesis** = stated, not tested. Scores are 1–5 (5 = best), each with its reason.

## 1. Central claim (one), and tonight's attempt to break it

**Claim.** An acoustic satellite of a Bragg reflection is a dynamical reflection with its own
amplitude. Written to first order in the strain, that amplitude gives a closed-form transfer function
H(ω; κ/2πα) between the strain spectrum and the time-Fourier weight of the reflectivity, valid
inside the Darwin width where the Sondhauss–Wark (2003) plane-wave form fails (Entries 59, 64, 65),
**for a single ray**: a beam whose spread in detuning exceeds the ringing width κ + π/x measures the
kinematic law times a flat factor κ√(π/2)/σ_η instead (Entry 82).
When the amplitude is not small the satellite forms dynamically with extinction length
Λ_P/|J_n(S)|, which is the poster's formation law (Entries 1, 9). The first-order amplitude ν is a
sufficient linearity criterion: ν ≲ 0.1 keeps the weight linear to 2 % (Entry 68) on the scanned
domain (Bragg, x ∈ {2, 5} outside the Darwin width, one packet family, weight at the spectral peak).

Where each part stands:

| part | status | evidence |
|---|---|---|
| first-order linearisation of TT for a phonon | prior art | Sondhauss & Wark 2003 eq. 25–29 (Entry 64) |
| kinematic heterodyne weighting W ∝ \|r0\| \|ε̂\|/ω | probably known in kind | Shayduk 2013 names the oscillation (Entry 60); not seen written as a transfer function |
| closed form H, Bloch-wave ansatz, valid inside the Darwin width | confirmed | TT to 1 % at κ/2πα = 0.12 incl. plateau; 0.5 % at 0.010 and 0.046 out to x = 100 (Entries 59, 64, 65) |
| S&W eq. 29 fails within ~1 Darwin width (factor 3–8) | confirmed, in our variables | Entry 64 §5 |
| H is a factor-2 bias for Si 400 at 17.5 keV over 20–100 GHz | confirmed **only for a long window** (the strain must travel ≥ max(1/κ, x/π) Λ during the window, whether it is a train or a ps packet followed for ≥ 50 ns); wrong for the 140 ps Chalise window, where H is flat to 4 % | Entries 65, 71, 75 |
| anomalous absorption χ_ih as a second crystal parameter | **rule fired** (Entry 77): with the X0h ratio χ_ih/χ_i0 = 0.93, H changes by 0.9 % (steady), 1.0 % (14 ns), 1.4 % (52 ns window) in the 20–150 GHz band at 17.48 keV and 3–4 % at 8 keV, opposite in sign on the two sides of the Bragg peak; the extended closed form (complex c) is TT-verified to 0.3 %; the ±η side asymmetry H(−η)/H(+η) = 1.02–1.09 is strain-spectrum-free and predicted, not verified | Entry 77 |
| ±η side asymmetry of H with χ_ih, strain-spectrum-free, reads χ_ih/χ_i0 | **confirmed in the model** (Entry 78, pre-registered): full TT vs closed form to 0.06 % at two energies, four windows; two packets agree to 0.1 %; size 1.9 % (17.48 keV, steady), 2.8 % (52 ns), 6.2 % (8 keV, 20 GHz); sensitivity dA/d(χ_ih/χ_i0) = 0.02–0.07; \|R0\| ratio 1.0000 in the phonon band (x ≥ 27), to be carried near the Darwin width | Entry 78 |
| asymmetric Bragg geometry | **confirmed, no solver change** (Entry 81): exact rescaling of the symmetric solver (TT to 1e-16), one geometry factor g → g(√b + 1/√b)/2 with b = γ_0/\|γ_h\|, window numbers scale with Λ_asym (T_90 52 → 33 ns at b = 4, 16.5 ns at b = 10 for Si 400 / 17.48 keV) | Entry 81 |
| beam divergence / bandwidth: the crystal factor is a single-ray property | **confirmed** (Entry 82): a beam spread in detuning σ ≥ 0.3 Darwin half-widths makes H_conv flat in f at every window, equal to κ√(π/2)/σ_η (TT beam average vs closed form 0.02 %); H(f), the window law, χ_ih and the ±η asymmetry need σ ≲ κ + π/x: 0.03–0.4 µrad and 2–4 meV at Si 400 | Entry 82 |
| beam condition at GHz-train frequencies (x = 1–10) | **confirmed, single-ray condition stands** (Entry 83): the 5 % tolerance allows σ ≤ 0.6 µrad / 7 meV (x = 1) to 0.3 µrad / 3.6 meV (x = 5) at 8 keV, ≤ 0.14 µrad / 9 meV at 17.48 keV; at x ≲ 2 the beam average first raises the factor 1.3–2× (resonance at π/x off the nominal detuning), so the flat law of Entry 82 holds only for x ≳ 10 | Entry 83 |
| beam condition for the formation law (P2, satellite centre) | **confirmed, condition table** (Entry 84): at σ_η = w_n (the satellite's own rocking width 2πα\|J_n\|) the first maximum loses 13–21 % (Bragg) and 25–40 % (Laue); both criteria (period 1 %, amplitude 5 %) hold only for σ_η ≤ 0.3 w_n in Bragg and 0.1–0.3 w_n in Laue, i.e. 0.1–0.3 µrad at Si 440 / 17.48 keV, 0.6–1.4 µrad at Si 400 / 8 keV; the integrated population is the beam-tolerant observable (its law not established) | Entry 84 |
| integrated satellite population as the beam-tolerant P2 observable | **rule failed as written** (Entry 86; the pre-registered amplitude metric compared differently normalised means, a stated flaw); post hoc, shape-normalised: the Laue population has the flat-beam shape from σ_η ≈ 3 w_n (1.3–3 µrad at Si 440 / 17.48 keV, one-channel detector 9–18 µrad), Bragg not even at 10 w_n; P_flat's shape is 13–38 % off sin²/∫J_0/tanh², so the population needs its own TT formation curve | Entry 86 |
| Laue integrated population as the beam-tolerant P2 observable (corrected metrics) | **confirmed with a condition** (Entry 87): rule at 3 w_n fails in 3 of 5 runs; all step runs pass from σ_η = 5 w_n (2–5 µrad at Si 440 / 17.48 keV, one-channel detector 9–18 µrad), the 1 t_tr ramp from 10 w_n (10 µrad); the law is the TT channel integral P_0(t) (tabulated), not sin² | Entry 87 |
| polarisation and static Debye–Waller | **absorbed** (Entry 88): π is the σ form with α → Cα, κ fixed (g/C, x/C, ν → Cν), exact; an unpolarised beam is the two-term closed-form sum, TT to 0.02 %; H_mix within 0.3 % of H_σ at Si 400 / 8 keV (C² = 0.12) but 20–32 % above at Si 440 / 17.48 keV (C = 0.73), so the polarisation state must be known there (vertical diffraction plane = σ); e^{−M} is already in the X0h constants | Entry 88 |
| closed form in the window length, H(f; T_w) | confirmed: the packet-kept first-order form (CF-B) matches TT to 0.3 % at every window 0.14–300 ns, two energies, three packet lengths; the packet-free short form (CF-A, Entry 71's hand estimate) is **rejected** below ~100 τ by the pre-registered 1 % rule | Entry 75 |
| signal maximum ≠ Fisher-information maximum (directive B2) | **negative** for a ps packet at one angle: amplitude maxima coincide (3 GHz apart), crystal moves nothing | Entry 70 |
| formation law I_n = sin²/tanh²(πα\|J_n\|t̃), no fit | confirmed | Entries 1, 9, 19, 24; tests |
| one dimensionless number joins the two regimes (directive B3) | **falsified in strong form** | Entry 68: spread 10–14× at fixed ν |
| ν ≲ 0.1 ⇒ linear to 2 % | confirmed as a bound on the scanned domain | Entry 68 §4 |
| Chalise et al. Fig. 2 tails = second-order strain terms | consistent in location; the low-frequency deficit is only partly explained (K = 0.65–0.88 vs their ~4×) once H is removed for their window | Entries 66, 71 |

**Falsification experiment run tonight** (`explore/nb/e68_nu_collapse.py`, Entry 68). Rule fixed
before the run: if 1 − K at equal ν spreads by more than ×2 across packet length, absorption and
satellite position, the "one criterion" form is false. Result: ×14 at ν = 0.1, ×10 at ν = 0.3.
So the claim is kept only in the form above (H is the correction, ν the guarantee); the sentence
"a single dimensionless parameter governs the crossover" must not appear in the paper.

**Second correction, same night (Entry 71).** H is the infinite-window limit. It needs the
measurement window, and so the strained depth, to exceed both the absorption depth 1/κ and the
extinction-phase depth x/π (tens of ns at Si 400): for the 140 ps packet of Chalise et al. the
crystal factor is flat to 4 % and the Entry 65 "factor 2" does not apply. H is therefore a property
of **long-window** measurements: the strain must travel ≥ max(1/κ, x/π) Λ during the window. A GHz
train read over tens of ns satisfies it, and so does a picosecond packet followed over a delay scan
of ≥ 50 ns (Entry 75: the 10 ps packet shows the full steady factor, H(20)/H(100) = 0.558 at
52 ns). The window condition therefore does not by itself tie P1 to the poster's apparatus
(physics cross-review, PR #16); sharing one apparatus remains an option, see below.

Reproduce:

```
uv run python explore/nb/e64_sw_compare.py            # S&W vs Bloch-wave first order vs TT (9 s)
uv run python explore/nb/e65_si400.py                  # H for Si 400 at 8 and 17.48 keV
uv run python explore/nb/e68_nu_collapse.py run stats fit2 fig   # tonight's falsification (28 s)
uv run python explore/nb/e70_fisher.py run fig                   # B2, one case (24 s)
uv run python explore/nb/e71_window_H.py                         # window condition for H (11 s)
uv run python explore/nb/e75_window_cf.py run posthoc fig        # H(f; T_w) closed form vs TT, boundaries (110 s)
uv run python explore/nb/e77_chi_ih.py run ttp ttm conv60 fig    # chi_ih: second parameter, side asymmetry (170 s)
uv run python explore/nb/e85_si440_case.py                       # Si 440 / 17.48 keV GHz-train case table (1 s)
```

**Fourth step, Day 2 (Entry 77).** The Bragg solver now takes χ_ih (`kappa_h`, f175d12). The
pre-registered 1 % rule fired at the 52 ns window (1.4 %), so H carries two tabulated crystal
parameters, g and χ_ih/χ_i0, no fit. New in kind: the two detuning sides of one phonon frequency
get different H (1.6–1.9 % at 17.48 keV, 3–9 % at 8 keV), a ratio in which the strain spectrum
cancels; the measured weight ratio also carries |R0(−η)|/|R0(+η)|, which is 1.0000 in the phonon
band (x ≥ 27) and 10 % near the Darwin width (physics cross-review). Pre-registered and confirmed in
the model as Entry 78 (T2f): TT and closed form agree to 0.06 %, two packets to 0.1 %, and the ratio
reads χ_ih/χ_i0 with sensitivity 0.02 (17.48 keV) to 0.07 (8 keV) per unit. This is the model's
experimental test of the second parameter, and a Fig. 3 panel candidate.

**Third step, Day 2 (Entry 75).** H(f; T_w) is now one closed form from the picosecond packet to the
steady train: at Si 400 / 17.48 keV the crystal factor is flat to 1 % up to 5 ns, departs by 10 % at
20 ns (T_flat), is halfway to the Entry 65 value at 33 ns and within 10 % of it at 52 ns (T_90); it
ripples around the steady value with period 2x Λ/v_L. The packet-free short form fails by 8 % at
140 ps (packet group delay and the −η image term) and was rejected under the rule; the packet-kept
form is what goes into the paper.

**Next falsification experiments, smallest first.** (T2a done: Entry 75. T2b done: Entry 76. T2c done: Entry 77.)
1. Anomalous absorption χ_ih in Bragg: done (Entry 77, `solve_bragg(kappa_h=)`); rule fired at the
   52 ns window (1.4 %), second parameter kept. Next smallest: pre-register the ±η side asymmetry
   as a claim and test it in full TT with a packet (not a single mode) at both energies. (T2c)
2. Asymmetric geometry (γ_0 ≠ γ_h): done, Entry 81 (T2g). No solver change: exact rescaling, one
   geometry factor in g, window condition in units of Λ_asym.
3. A smooth strain front instead of the step front of the Entry 58 packet: changes the high-frequency
   wing, and may change the saturation route of Entry 68.
4. Read Chukhovskii & Förster 1995 and Wark & Lee 1999 (time-dependent TT) for an earlier Bloch-wave
   first-order form. If one exists, the H part drops to "application"; the ν bound and the Si 400
   number stand. (T2d)
5. B1 degeneracy at larger S (Entry 72 iii): done, class holds to S = 6 (Entry 76, T2b).

## 2. Candidates, scored

Scores: N = academic novelty, E = experimental importance, D = implementation difficulty (5 = easy),
G = chance of completion before graduation (2027 spring).

| # | candidate | N | E | D | G | one-line reason |
|---|---|---|---|---|---|---|
| P1 | Transfer function H(f; T_w; g, χ_ih/χ_i0) + window condition + linearity bound for reading phonon spectra from time-resolved Bragg reflectivity (Entries 59, 64–66, 68, 71, 75, 77, 82) | 3 | 2 | 4 | 5 | Closed form done and TT-verified; after Entry 82 the f-dependence of H is a single-ray property seen only with meV / sub-µrad resolution, while any ordinary beam measures the kinematic law times a flat κ√(π/2)/σ_η (E lowered 3 → 2); the practical message for ps ultrasonics is "no crystal factor, divide by |R0|/ω"; novelty bounded by S&W 2003 and two unread 1990s papers |
| P2 | Temporal Pendellösung of satellites behind a travelling front: formation law, Bessel ladder, period-ratio calibration of S (Entries 1, 9, 13, 19, 24) | 3 | 3 | 5 | 4 | Clean, no-fit, both geometries, tests; prior art DeCamp 2001 (n = 0 pulse) + Polikarpov 1987/88 (static extinction length) makes it "extension to every order"; needs a bulk GHz shear-wave experiment that does not yet exist; the elastic side is clean for all S ≤ 6 at 1–2 GHz (Entry 80), so the only open feasibility item is the transducer peak power, which keeps E at 3 until a number is in hand |
| P3 | Path memory: Θ_n law, transient-only satellite, overshoot map (Entries 28, 30, 34, 67) | 2 | 3 | 5 | 4 | Robust in the model, but at α = 0.5 it is kinematic diffraction from a graded grating and the slow-ramp limit is the Bessel ladder (Entry 34); best as a section of P2, not a paper |
| P4 | Overlap regime: Bessel-zero resurrection with a computable third-order coupling and shifted zero (Entries 45, 51, 52, 54) | 4 | 1 | 3 | 2 | Possibly new, but observability is 4e-5 of the window integral (Entry 54) and Entin / Iolin–Entin are unread; theory only |
| P5 | Second-order strain terms in the Chalise spectrum (Entry 66) | 2 | 3 | 5 | 5 | Known in kind (Bessel expansion); its value is the size at their strain and the fluence-series criterion; a section of P1 |

**What carries the paper now (after Entries 71, 75, 82).** Not the Chalise correction. It is the pair
"closed-form H valid inside the Darwin width + its window condition" together with the formation law,
**with the Entry 82 caveat**: the shape of H and its window law belong to a single ray; a beam whose
spread in detuning exceeds the ringing width κ + π/x (0.03 µrad / 2 meV at 17.48 keV, 0.4 µrad / 4 meV
at 8 keV) measures a flat crystal factor κ√(π/2)/σ_η at every window, i.e. the kinematic law. The
long-window H is therefore an inelastic-scattering-grade measurement, and the ps-ultrasonics message
is the simpler one (no crystal factor for an ordinary beam).
H needs a long window: the strain travels ≥ max(1/κ, x/π) Λ during it, which for Si 400 at 17.48 keV
is ≥ 22 ns (strain travels ≥ 185 µm during the window). Two measurements satisfy it: (i) a
picosecond-ultrasonics packet with a delay scan of ≥ 50 ns, i.e. P1 has an experiment of its own
that needs no GHz train; (ii) a GHz bulk train read over tens of ns, the poster's Si 440 Laue
geometry (t_tr 28–51 ns, Entry 20), where the formation law of P2 appears in the same data. The
"one apparatus" of the recommendation is option (ii), a choice, not a consequence of the window
condition (physics cross-review). Entry 75 gives the window numbers for Si 400 / 17.48 keV, 20 GHz against 100 GHz:
the crystal factor is flat (within 10 %) up to **T_flat = 19.6 ns**, halfway to the steady Entry 65
value at **T_half = 32.8 ns**, and within 10 % of it from **T_90 = 52.2 ns** on (TT-confirmed to
0.002). A measurement with a 20–50 ns window sits on the steep part, where H depends on the window
length itself; the experiment should therefore read either ≤ 5 ns (H flat to 1 %, kinematic
estimator) or ≥ 52 ns (steady H), and in between must know its window to 10 %. Without such a measurement the paper is theory with a proposal. T1c (Entry 73): on
the elastic and TT side S ≤ 6 at 1–2 GHz is clean (strain ≤ 2e-4; the second-harmonic column of
Entry 73 was withdrawn in Entry 80: a pure shear wave along ⟨100⟩ has no collinear second harmonic,
and the non-phase-matched longitudinal 2ω component is ≤ 1e-3 and invisible to the in-plane
reflection); the drive is 6–100 W/cm² peak acoustic intensity for S = 1–2 and 230–900 W/cm² for
S = 3–6, pulsed, and no verified GHz bulk-transducer number is in hand. A real
transducer rings up over 0.2–4 t_tr, so the observable regime is the ramp map (Entry 67), not the
step collapse.

**Recommendation (revised after Entry 82).** The theory paper P1 + P2 stands, but its experimental
section must lead with the beam condition: for σ ≫ κ the estimator needs no crystal factor and the
window law is unobservable; the single-ray H, χ_ih and the ±η asymmetry are proposed as a meV /
sub-µrad measurement (or an angle-resolving detector plus a meV monochromator). Previous text, kept
for the single-ray case: one paper, P1 + P2, around one experiment: a GHz bulk shear wave entering a
Si crystal, read over tens of ns. Its linear regime (ν ≲ 0.1) is the H-corrected spectrometer, its
dynamical regime is the formation law, and a long window lets both appear in the same measurement.
The alternative, P1 alone on a ps packet with a ≥ 50 ns delay scan (an existing kind of measurement,
Reis 2001 and the Chalise geometry), is open and closer to current practice; it leaves P2 as the
poster's own paper. P5 is a subsection on picosecond packets (where only K matters).
P3 stays inside P2. P4 is a future project (Entry 51 F). Before tonight the recommendation was P1
alone on the strength of the Chalise correction; Entry 71 removed that strength.

Risks, in order:
- Entry 82: the long-window crystal factor needs σ ≲ 0.03–0.4 µrad and ≲ 2–4 meV; with an ordinary
  beam the measurable crystal factor is flat. This is now the first risk. Entry 83: the GHz-train
  band relaxes this only to 0.3–0.6 µrad / 4–7 meV at 8 keV, so the E scores of P2 and P1+P2 are
  not raised. Entry 84: the P2 centre observable needs σ_η ≲ 0.3 w_n (0.1–0.3 µrad in the poster's
  Laue geometry), so P2 is a sub-µrad measurement too as a centre observable. Entry 86: the Laue
  integrated population is beam-tolerant from σ_η ≈ 3 w_n (1.3–3 µrad, one-channel detector), which
  would lift the condition for the Laue train, but only under a corrected pre-registration (the one
  run had a mis-specified amplitude metric) and with a TT-tabulated formation curve (13–25 % off
  ∫J_0); Bragg populations are not clean. Entry 87 (corrected metrics): the condition is σ_η ≈ 5 w_n
  (2–5 µrad) for a step front and 10 w_n (10 µrad) for a 1 t_tr ramp, an ordinary high-resolution
  beam; E of P2 can go to 4 if Isaac accepts that condition (the pre-registered 3 w_n failed, so it
  is not raised here).
- No long-window (tens of ns) time-resolved satellite measurement exists yet; the paper's
  experimental relevance is a proposal. Entries 73, 80: elastic side clean for all S ≤ 6; the
  transducer peak power (6–100 W/cm² at 1–2 GHz for S = 1–2, up to 900 W/cm² at S = 6) is unverified.
- A referee calls the Bloch-wave correction of S&W routine. Answer: the paper rests on the window
  condition, the bound and the formation law, not on the linearisation.
- An earlier first-order Bloch-wave form exists (Chukhovskii & Förster, Wark & Lee). Not read.
- Papers a referee will name (physics cross-review, bibliography unchecked, on the T2d list): for
  P1 Reis et al. PRL 86, 3072 (2001), Trigo et al. Nat. Phys. 9, 790 (2013), Schick et al.
  udkm1Dsim CPC 185, 651 (2014); for P2 Bucksbaum & Merlin SSC (1999), Zolotoyabko and Iolin on
  acoustically modulated crystals, Grigoriev et al. APL (2006). The novelty sentences are to be
  written against these.
- Chalise's energy and normalisation are unknown; the comparison stays "same sign, same order".
- Model limits: σ or π by rescaling, unpolarised as a two-term sum (Entry 88), normal + anomalous absorption in Bragg (Entry 77; Laue still lossless), symmetric or asymmetric Bragg by rescaling (Entry 81; Laue symmetric only), plane wave; beam divergence and bandwidth handled as a single-ray average (Entries 82–87); lateral acoustic uniformity assumed.

## 3. Figure plan (figures go to `figures/paper/`; poster figures are frozen)

1. **Fig. 1 — the estimator and the crystal factor.** (a) R(t) at one detuning for a packet, with its
   Fourier component at ω = |η| marked; (b) W/W_kin against ω/ω_D: the one-curve collapse for three
   packet spectra at κ/2πα = 0.12 (`e59_transfer.py fig`, `figures/chalise_transfer.png` b).
2. **Fig. 2 — where the plane-wave first order fails.** |F_pred|/|F_TT| at the satellite against
   ω/ω_D for S&W eq. 29, S&W normalised, and the Bloch-wave (Riccati) form; the near-main-peak
   error alongside (`e64_sw.json`). Message: the Bloch-wave depth dependence is what keeps first
   order valid inside the Darwin width.
3. **Fig. 3 — a real crystal, and the window condition.** H(f; T_w) against T_w for Si 400 at
   17.48 keV, closed form (lines) with TT (circles), and the ratio H(20)/H(100) against T_w with
   T_flat, T_half, T_90 marked (`figures/paper/e75_window_H.png`, Entry 75): the suppression exists
   only when the window exceeds 1/κ and x/π. One panel with K(f) for
   a picosecond packet (`figures/si400_harmonic.png`), where only the second-order term is left.
4. **Fig. 4 — the linearity bound and the dynamical satellite.** (a) 1 − K against ν for the 24
   series of Entry 68 with the tanh² and sin² laws (`figures/paper/e68_nu_collapse.png`, left panel);
   (b) the formation law of the poster, I_n(t̃) = sin²(πα|J_n|t̃) collapse (`figures/poster/f1_collapse.png`
   re-drawn in `figures/paper/`), as the ν ~ 1 limit.

## 3b. Directive questions, where each stands (T1b)

| question | entry | result |
|---|---|---|
| B1 two fields, same steady state, different transient | 72, 76 | confirmed, and bounded: order n reads the t_tr running mean of J_n(S A); zero-mean perturbations are degenerate at that order (TT 9e-7 at S = 1.5; ≤ 3e-6 for S ≤ 6, both orders, both J_1 branches, Entry 76); the other order separates the same pairs by 1e-3 to 1.4e-2 (law and TT agree to 9e-6); the class erodes with the coupling ratio R (2.9e-4 at R = 0.1), not with S |
| B2 signal maximum vs Fisher maximum | 70 | negative for a ps packet at one angle (3 GHz apart, crystal moves nothing); open for long trains |
| B3 one dimensionless criterion | 68, 71 | falsified in strong form; ν ≲ 0.1 bound and the window condition survive |
| B4 features when independent satellites fail | 45, 51, 52, 54 (no rerun) | Bessel-zero resurrection, confirmed in model, unobservable (4e-5), "possibly new" pending Entin / Iolin–Entin |

## 4. Open decisions for Isaac

- Paper = P1 (transfer function) with P2 inside, or P2 (formation) as the first paper? P1 is
  closer to finished and tied to a current measurement; P2 is the poster's identity.
- χ_ih is in `solve_bragg` (Entry 77, PR #16). Asymmetric geometry needs no solver change (Entry 81):
  closed. Whether an asymmetric cut (b = 4–10, T_90 = 33–16 ns) goes into the proposal as the way
  to shorten the delay scan.
- Whether the ±η side asymmetry (Entry 78, model-confirmed, unobserved) goes into the paper as a
  proposed measurement of χ_ih/χ_i0.
- Library access for Chukhovskii & Förster 1995 and Wark & Lee 1999.
- Whether to tell Chalise et al. that, in this model, their low-frequency deficit is not the crystal
  factor (Entry 71) and only partly the second-order term (Entry 66); and whether to correct
  `notes/qa_prep.md` §6b (the poster itself is frozen).
