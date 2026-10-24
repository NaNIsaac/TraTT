# Referee pass on the draft skeleton (BACKLOG T2j, draft tier; no new computation)

Written 2026-10-23 against `notes/paper_draft.md` (d20df9c). Three objections a referee is most
likely to raise, in referee wording, each followed by what already answers it (entry, number) and
what does not yet exist: the smallest computation or measurement that would. The three are the
three axes boss named: prior art, model limits, experimental feasibility. A fourth, shorter list
of secondary objections follows.

## Objection 1 (prior art). "The first-order treatment of a phonon in dynamical diffraction is in Sondhauss & Wark (2003), and time-dependent Takagi–Taupin codes (udkm1Dsim; Chukhovskii & Förster; Wark & Lee) have integrated strain profiles through the Darwin width for twenty years. What is new beyond rewriting the known linear response with a window?"

Referee's strongest form: the Bloch-wave depth dependence that keeps first order valid inside the
Darwin width is a routine consequence of the Riccati form; the window dependence is a routine
finite-time Fourier integral; neither is a physical result.

**What answers it now.**
- Entry 69 aligns S&W eq. 25–29 term by term with Entry 59: S&W's form is the plane-wave
  (constant-amplitude) ansatz and is off by a factor 3–8 within one Darwin width (Fig. 2,
  `e64_sw.json`); the Bloch-wave form is exact to first order at every detuning (1 % at
  κ/2πα = 0.12, 0.5 % at 0.010 and 0.046 to x = 100, Entry 65 §3 converged). The novelty
  sentence can be written as "the plane-wave first order fails where the satellite sits for
  low-order Si reflections (x ≲ 5), and the Bloch-wave first order does not"; that is a statement
  about where a published formula can be used, not a rewriting.
- The window condition is not a Fourier-integral triviality: Entry 75 shows the *packet-free*
  closed form (the naive truncated kernel) fails by 8 % at short windows and only the packet-kept
  form holds; the design numbers T_flat/T_half/T_90 and the e^{−κT} ripple (TT-confirmed) are
  the physical content, and Entry 71 shows a published-type measurement (140 ps Chalise window)
  sees none of the steady factor that Entry 65 would have applied to it.
- Entries 77–78 add a second parameter and a measurable consequence (±η asymmetry) that no
  linear-response paper we have read states.

**What is missing, smallest first.**
1. Reading Chukhovskii & Förster 1995 and Wark & Lee 1999 in full (T2d, PDFs pending). If either
   writes the first-order Bloch-wave amplitude in closed form, claim 2 of the draft drops to
   "application" and the paper leans on claims 3, 4, 6, 7, 12 (window, design numbers, χ_ih,
   asymmetry, degeneracy), which are not in that literature as far as read.
2. One paragraph placing H against udkm1Dsim: a numerical TT code gives H by brute force for one
   case; the paper gives it in closed form with its two parameters and its window law. State
   that the closed form was checked against our own TT, not against udkm1Dsim; a cross-check
   against udkm1Dsim on one case (Si 400, 17.48 keV, the Entry 66 packet) is the minimal
   computation that would pre-empt the "your TT is wrong" variant of this objection. Not run.
3. Reis et al. 2001 and Trigo et al. 2013: cite as the two experimental settings (ps packet with
   dynamical interpretation; Fourier-domain phonon spectroscopy) and say explicitly that neither
   corrects by a window-dependent crystal factor. Reading them is a T2d item.

## Objection 2 (model limits). "A plane-wave, monochromatic, σ-polarised, two-beam, quasi-static model with the strain entering only as a phase. Real measurements have a divergent, polychromatic beam, a finite footprint, both polarisations, thermal Debye–Waller change under GHz drive, and a transducer-launched wave that is neither a step nor a perfect plane. Which of the paper's percent-level numbers survive?"

Referee's strongest form: the 0.9–1.4 % χ_ih effect, the 1.6–6 % ±η asymmetry and the 10 %
window-knowledge requirement are all smaller than or comparable to what a beam of 10 µrad
divergence integrates over the Darwin width.

**What answers it now.**
- Divergence and bandwidth enter as an average of H over detuning. In the phonon band (x ≥ 27)
  H varies slowly with x (Entry 65 table: 0.786 → 0.909 from 40 to 60 GHz), so an angular spread
  of a few Darwin widths changes H by a few per cent and is itself a known instrument function;
  the paper must say H is to be convolved, and the convolution does not move the window numbers
  (they come from κ and x/π, Entry 75 §1), nor the sign of the ±η asymmetry (it is odd in η,
  Entry 78, so any symmetric angular average preserves it).
- π polarisation: done, Entry 88. The π component is the σ form with α → Cα and κ fixed (exact);
  an unpolarised beam is the two-term closed-form sum, verified against TT to 0.02 %. At Si 400 /
  8 keV (C² = 0.12) H_mix is within 0.3 % of H_σ; at Si 440 / 17.48 keV (C = 0.73) it is 20–32 %
  above, so the polarisation state must be stated there (vertical diffraction plane = σ).
- Quasi-static: β_qs ≤ 1.2e-2 at 1–2 GHz (Entry 73); for 100 GHz packets the retardation is
  still 1e-2 over 10 Λ (Entry 62). Covered.
- Strain only as a phase: DESIGN.md item 7; the elasto-optic change of χ_h at 1e-5 strain is
  of order the strain itself (1e-5 relative), below every number quoted.
- Step vs transducer ramp: Entry 73 §3 and Entry 67 map the ramp regime; the step is not claimed
  as the observable regime.
- The ±η asymmetry and the χ_ih correction are ratios in which the strain spectrum cancels
  (Entry 78: two packets agree to 0.1 %), so they are insensitive to the packet shape, which is
  the part of the model a referee distrusts most.

**What is missing, smallest first.**
1. (done, Entry 82, T2k) The convolution was run, in full TT and closed form, and the referee is
   right in a stronger way than "a few per cent": the ringing linewidth is κ + π/x, 0.01–0.05 Darwin
   half-widths, so any beam spread of 0.3 half-widths or more dephases the ringing after 1/σ_η and
   makes the crystal factor flat, κ√(π/2)/σ_η, at every window. The answer in the paper: (a) for an
   ordinary beam the estimator needs no crystal factor, which is the practical result; (b) H(f; T_w),
   χ_ih and the ±η asymmetry are single-ray observables proposed for a meV / sub-µrad measurement.
   The "slowly varying in the band" argument above is withdrawn: it assumed the ray's own H, not the
   beam-averaged ringing.
2. Debye–Waller: done, Entry 88. The static e^{−M} (0.94 Si 400, 0.88 Si 440 at 300 K) multiplies
   χ_rh and χ_ih and is already inside the X0h constants used throughout; S = h·u is untouched and
   the coherent displacement enters through J_n(S). Drive heating at the 1e-4 duty cycle (mW/cm²)
   does not move M at any quoted precision. A sentence in the limitations.
3. Two-beam: multiple-beam cases near Si 400 at 17.48 keV should be listed as excluded azimuths.
   No computation; a sentence.

## Objection 3 (experimental feasibility). "No measurement of the kind proposed exists. The long-window regime needs the strain to travel 185–440 µm during the window, a GHz bulk transducer at 6–900 W/cm² peak, and a window known to 10 %. The ps-packet route needs a 50 ns delay scan at 1 % reflectivity precision. Either way the paper is a theory with a proposal, and the proposal's key number (transducer power) is unverified."

**What answers it now.**
- Entries 73 and 80: the elastic side is clean for all S ≤ 6 at 1–2 GHz (strain ≤ 2e-4, no
  shear second harmonic, longitudinal 2ω component ≤ 1e-3 and invisible to the in-plane
  reflection); drive power is the only limit, and S ≤ 2 needs 6–100 W/cm² peak at a 1e-4 duty
  cycle (mW/cm² average).
- Entry 75: the ps-packet route does not need the train; a 10 ps packet followed to 52 ns shows
  the full factor (H(20)/H(100) = 0.558). Existing ps-ultrasonics X-ray setups (Reis-type) scan
  delays of that length; what is new is only the required length and the 1 % precision on the
  ratio of two frequencies' Fourier weights.
- Entry 81: an asymmetric cut shortens T_90 to 33 ns (b = 4) or 16.5 ns (b = 10) at the price of
  the γ_0/|γ_h| intensity factor; this is the knob that brings the steady regime into a shorter scan.
- Entry 78: the ±η asymmetry is a *ratio* measurement at one delay, free of the strain spectrum and
  of |R0| in the phonon band, so it does not need absolute calibration; its size (2–6 %) is above
  the 1 % precision asked above.
- Entry 20 (X0h constants) and Entry 73 give every X-ray number; nothing there is estimated.

**What is missing, smallest first.**
1. A transducer number: peak acoustic intensity of a GHz bulk shear transducer on Si (ZnO or AlN
   film, or a bonded LiNbO3 plate). Library/Isaac item (DotFiles#57); without it the train route
   stays "feasible on the elastic side, unverified on the drive".
2. A signal-to-noise estimate for the ratio measurement: with R0 ≈ πα/η in the tail, the
   oscillation amplitude is 2|R0||δ| ≈ 2(πα/η)² S|G|... at 1e-5 strain (Fig. 1a: 1 % of |R0|² at
   40 GHz), so the ±η ratio to 0.1 % needs 1e6 photons per detuning per delay bin; a one-line
   count against a synchrotron/XFEL flux belongs in §5.5. Closed-form arithmetic, not run.
3. The 10 % window knowledge: in a delay scan the window is the scan length, known exactly; in a
   gated detector it is the gate. State which, and that T_flat/T_90 move with x and g by Entry 75.
4. Ideally one measurement. The smallest: a ps-packet delay scan on Si 400 at 8 keV (1/κ = 4.8 ns,
   steady from 30 ns, H(20)/H(100) = 0.93 and the ±η asymmetry 6 % at 20 GHz, Entries 71, 78),
   where the effects are largest and the window shortest. This is the experiment to propose first.

## Secondary objections (one line each)

- "ν ≲ 0.1 is only shown for x ∈ {2, 5}, one packet family." Yes; the bound is stated on that
  domain (Entry 68); extending it to x = 30–200 (the Si 400 band) is a 30 s scan not run.
- "The degeneracy class is a property of your approximate law, not of TT." It is checked in full
  TT to 3e-6 (Entries 72, 76) and its R boundary found at 0.2 (Entry 79); the law's own error is
  ten times larger, so the class is a TT property.
- "χ_ih/χ_i0 = 0.93 is assumed equal at both energies." X0h gives 0.925 and 0.931 (Entry 20
  script); use each.
- "The Chalise comparison is unquantified." Agreed and stated (Entries 57, 66, 71); the paper
  says "same sign, same order" and that their 140 ps window sees no crystal factor.
- "Bragg only; the poster is Laue." The formation law is both geometries (Fig. 4b); H and χ_ih are
  Bragg; Laue absorption is excluded (DESIGN.md item 4) and said so.
- "Why no instrument function at all?" Objection 2, item 1 is the minimal fix.

## Summary for the authors

| objection | strongest answer in hand | smallest missing item |
|---|---|---|
| prior art | Entries 64, 69 (S&W fails at x ≲ 5), 75 (window law is not trivial), 77–78 (new parameter and observable) | read C&F 1995, W&L 1999 (T2d); one udkm1Dsim cross-check |
| model limits | Entries 82–88: beam spread → single-ray average (flat factor at x ≫ 1, non-monotonic at GHz); centre law needs 0.3 w_n, Laue integrated population 5 w_n (2–5 µrad); polarisation and Debye–Waller absorbed as rescalings | lateral acoustic uniformity (transducer width, not modelled) |
| feasibility | Entries 73/80 (elastic clean), 75 (ps route), 81 (asymmetric cut), 78 (ratio at one delay) | transducer number (DotFiles#57); photon-count line; propose Si 400 at 8 keV first |
