# First proposed experiment, one page (BACKLOG T2l, draft; posting is Isaac's; no new computation)

Written 2026-10-23 from Entries 65, 71, 73, 75, 77, 78, 80, 81, 82. Numbers are the model's
(symmetric Bragg, σ polarisation, plane-wave rays, first order in strain); nothing here has been
measured. The experiment: **Si 400 at 8 keV, a picosecond strain packet from the surface, time-
resolved reflectivity over a delay scan of ≥ 30 ns, at two detunings ±ω of the Bragg peak.**

## 0. The beam condition comes first (Entry 82)

Each incident ray rings, after the packet, at its own detuning frequency with damping κ; the
Fourier weight at a nominal phonon frequency collects only rays within that linewidth. At Si 400
/ 8 keV the linewidth is 0.046 Darwin half-widths = **0.38 µrad in angle, 4.4 meV in energy**
(Darwin half-width 8.3 µrad). Two regimes, and the proposal must say which it is in:

| regime | beam spread in detuning | what is measured | what it tests |
|---|---|---|---|
| A, ordinary beam | ≳ 0.3 Darwin half-widths (≳ 2.5 µrad or ≳ 25 meV) | kinematic law × a flat factor κ√(π/2)/σ_η, at every window; no crystal factor, no window law, no ±ω asymmetry | the estimator "divide by \|R0\|/ω" and the flat factor's scaling with the beam's σ |
| B, single-ray | ≲ 0.4 µrad **and** ≲ 4 meV (angle-resolving detector + meV monochromator) | H(f; T_w), its window law, the χ_ih correction, the ±ω asymmetry | Entries 65, 75, 77, 78 |

σ_η is the beam's spread in the deviation parameter η ∝ Δθ − (ΔE/E) tan θ_B, i.e. the projection of
the joint (θ, E) distribution, not the quadrature sum of divergence and bandwidth (physics
cross-review): with a Si 111 monochromator and a Si 400 sample the arrangement is dispersive and
σ_η is large; a non-dispersive (+,−) arrangement with matched d-spacing cancels the bandwidth in η
and leaves the first crystal's intrinsic width; an asymmetric (+,−) collimator with |b| ≪ 1 gives
σ_η ≈ √|b| Darwin half-widths (b = 1/50 → 0.14) with routine optics, which reaches the single-ray
condition only for strongly absorbing reflections (g ≳ 0.1). σ_η must be computed from the optics
for the chosen beamline before the regime is assigned.

Regime A is a routine measurement and gives the practical result for picosecond ultrasonics.
Regime B is an inelastic-scattering-grade beam (meV monochromator, 1e9–1e10 ph/s) with a 0.4 µrad
collimation or a detector that resolves the exit angle to that level; it is the one that measures
the paper's crystal factor. Where such beams exist (to be checked before writing; confidence in
brackets): meV monochromators at IXS beamlines, NSLS-II 10-ID at 9.1 keV (~1 meV) [medium],
ESRF ID28 / APS 30-ID HERIX / SPring-8 BL35XU at 17.8–23.7 keV (1.5–3 meV) [high, but these are
storage rings with 50–100 ps bunches, too long for the 20–150 GHz band unless a streak camera
(~1 ps) is the detector]; at an XFEL a four-bounce high-resolution monochromator reaches meV at
8–9 keV at ~1e-4 of the SASE bandwidth, i.e. 1e8 ph/pulse, 1e10 ph/s at 120 Hz [medium], which is
the flux assumed below; the 0.4 µrad angular condition is the harder one (XFEL divergence 1–2 µrad):
either a collimating optic or exit-angle resolution of 0.4 µrad (40 µm pixels at 100 m) [open].
The combination "meV + ps" is uncommon; it is the main risk of regime B.

A third option was checked (Entry 83): at GHz-train frequencies (x = 1–5) the ringing half-width
grows to g + 1/(2x) Darwin half-widths, but the 5 % tolerance on the single-ray factor still allows
only σ ≤ 0.6 µrad / 7 meV at x = 1 and 0.4 µrad / 5 meV at x = 2 (8 keV; 0.14 µrad / 9 meV at x = 2
for 17.48 keV), because the single-ray resonance sits π/x away from the nominal detuning and a beam of
a few tenths of a half-width lands on it (the factor first rises 1.3–2×, then falls). So the GHz-train
apparatus eases the energy condition slightly and the angular one not at all; it does not reach an
ordinary beam for the *single-ray* H. What an ordinary beam does see at GHz frequencies (x ≲ 10) is
a frequency-dependent, non-monotonic beam-averaged factor H_conv(x; σ_η) = p(π/x)/p(0) × (single-ray
shape), computable by the CF-B convolution (TT to 0.24 %, Entries 82–83) once σ_η is known from the
optics; that factor is itself a falsifiable prediction and the observable P1 + P2 share in the GHz
regime (pre-registration candidate, not yet run as such). For the formation law itself (satellite centre at S ~ 1, P2) the divergence question is
the satellite's own rocking width and is not computed.

## 1. Beamline requirements

- Time structure: an X-ray probe short against the phonon period (≤ 2 ps for 100 GHz; XFEL or
  streak camera; the 2–4 meV of regime B is transform-compatible with pulses ≥ 0.9 ps, and
  synchrotron bunches of 50–100 ps cannot resolve 20 GHz in any regime), delay range ≥ 30 ns in steps ≤ 2 ps near the packet and ≤ 50 ps later (the
  Fourier weight at 20–150 GHz needs the fine steps only inside the ringing; a two-step scan).
- Pump: a ps-ultrasonics transducer film giving a strain packet of peak strain ~1e-5 and carrier
  20–150 GHz (the Entry 66 packet: 87 GHz, 0.9 cycles); the linear regime ν ≲ 0.1 is then
  guaranteed (Entry 68), and the strain spectrum is read from a thin-film reference as usual.
- Angular/energy resolution: regime A needs nothing special; regime B needs ≲ 0.4 µrad and
  ≲ 4 meV at the sample (see 0).
- Detuning: two angles at ±ω/(dη/dθ) from the refraction-corrected Bragg angle, for ω = 20–60 GHz
  (x = 38–115 Darwin half-widths at 8 keV): in angle ±(ω/f_D) × 8.3 µrad with f_D = 0.52 GHz,
  i.e. ±0.3–1 mrad; the static rocking curve at the same angles is recorded for |R0|.
- Crystal: thickness ≥ 0.3 mm (the packet is 0.25 mm deep at 30 ns), (001) Si, polished.
- Flux (regime B, arithmetic from Fig. 1a scaled): the oscillation is 1 % of |R0|² at 1e-5 strain,
  |R0|² ≈ 1/(4x²) = 1.7e-4 at 20 GHz; a 1 % measurement of the Fourier weight needs ~1e4 relative
  precision per delay bin, i.e. 1e8 detected photons per bin, 100 bins, two detunings: 2e10
  detected, 1e14 incident at 1.7e-4 reflectivity: 3 h at 1e10 ph/s. The ±ω ratio (6 %) needs
  0.5 % and so 100× fewer. Regime A at 1e12 ph/s: minutes.

## 2. Expected signal (model numbers, Si 400 / 8 keV unless stated)

- Regime A: the Fourier weight at ω over the window equals |G(ω)| |R0(ω)| πS/κ × κ√(π/2)/σ_η,
  flat in ω to 0.3 % (Entry 82 at σ = 0.3 half-widths); the factor scales as 1/σ_η (vary a slit
  or the monochromator bandwidth by 3× and the weight changes 3×).
- Regime B, window law (Entries 71, 75): H(20)/H(100) = 1.014 at 140 ps, 0.998 at 1.4 ns, 0.958
  at 14 ns, 0.929 steady (from 30 ns on); the absolute H rises from 0.03 at 140 ps to 0.93 at
  30 ns (17.48 keV table of Entry 75 scaled by 1/κ; the 8 keV window table was not run, only the
  three points of Entry 71).
- Regime B, χ_ih (Entry 77): the steady H changes by −2.8 % (η > 0) and +3.2 % (η < 0) at 20 GHz
  with χ_ih/χ_i0 = 0.925; 4.5 % at 10 GHz.
- Regime B, ±ω asymmetry (Entry 78): H(−ω)/H(+ω) − 1 = 6.2 % at 20 GHz, 3.4 % at 40, 1.4 % at 100
  (steady; 5.8 % / 3.0 % / 1.2 % at 14 ns; 1.0 % at 1.4 ns), independent of the packet spectrum
  to 0.1 %, |R0| ratio 1.0000 in this band; sensitivity dA/d(χ_ih/χ_i0) = 0.068 at 20 GHz, so a
  0.1 % ratio resolves χ_ih/χ_i0 to 0.015.
- Asymmetric cut (Entry 81, computed at 17.48 keV, same formula at 8 keV): b = γ_0/|γ_h| = 4
  shortens the steady onset by 37 % and b = 10 by 68 %, at the cost of a 2.7° exit angle and the
  γ_0/|γ_h| intensity factor; an option if the delay range is the limit, not needed at 8 keV.

## 3. Three outcomes that would show the model is wrong

1. **Regime A shows a frequency-dependent crystal factor.** If, with an ordinary beam
   (σ ≥ 2.5 µrad or ≥ 25 meV), the Fourier-weight spectrum divided by |G(ω)| |R0(ω)|/ω is not flat
   to a few per cent across 20–100 GHz, or depends on the window beyond the coherence time
   1/σ_η, Entry 82's dephasing picture is wrong.
2. **Regime B shows no window law.** If with ≤ 0.4 µrad and ≤ 4 meV the ratio H(20)/H(100) stays
   at 1 ± 0.02 from 1 ns to 30 ns, the single-ray crystal factor of Entries 65/75 (and the ringing
   at the detuning frequency) is wrong; if it reaches 0.93 but at a window 2× different from
   30 ns, the absorption/extinction depths in the model are off.
3. **Regime B is symmetric in ±ω.** If H(−ω)/H(+ω) = 1 ± 0.01 at 20 GHz beyond 14 ns, the
   anomalous-absorption term of Entry 77 is wrong or χ_ih/χ_i0 is not 0.93; if the asymmetry
   has the opposite sign, the sign convention of Entry 77 (higher static reflectivity on the
   low-absorption edge) must be checked against the static rocking curve, which is recorded.

Any of the three is reported as a falsification, not fitted away (BACKLOG rule 20).

## 4. Transducer alternative (the poster's regime)

A GHz bulk shear train in Si 440 Laue or Si 400 Bragg, read over tens of ns (Entries 20, 73, 80):
S ≤ 2 needs 6–100 W/cm² peak acoustic intensity at 1–2 GHz, pulsed at 1e-4 duty cycle; the
elastic side is clean for all S ≤ 6 (no shear second harmonic); the only unverified number is the
transducer's peak output (DotFiles#57). It measures the formation law (P2) at S ~ 1–6, a different
observable (satellite centre intensity) whose divergence sensitivity is the satellite's own
rocking width w_n = 2πα|J_n(S)|: Entry 84 finds the centre law intact only for σ_η ≲ 0.3 w_n,
0.1–0.3 µrad at Si 440 / 17.48 keV (0.6–1.4 µrad at Si 400 / 8 keV Bragg, the more tolerant
geometry), so as a centre observable it needs the same sub-µrad beam as regime B; with an ordinary
beam the measurable quantity is the integrated satellite population: in Laue its shape is the
flat-beam one from σ_η ≈ 5 w_n (2–5 µrad at Si 440 / 17.48 keV; 10 µrad for a 1 t_tr ramp) with a
detector accepting one satellite channel (9–18 µrad), its formation curve is the TT channel integral
(tabulated in `e87_population_laue.json`; 13–25 % off the ∫J_0 law, so not a closed form), and in
Bragg it is not clean at any beam (Entries 86–87). One case in numbers (Entry 85, Si 440 / 17.48 keV, L = 163 µm, t_tr = 28 ns): the
linear regime of the filled crystal needs S ≤ 0.02 (strain ≤ 7e-7, mW/cm²), the formation regime
S = 0.5–2 (strain 1–7e-5, 1.6–100 W/cm², τ_1 = 8–19 ns); the Bragg window numbers of the same
reflection (longitudinal [110] wave, since a shear wave along the normal does not couple in
symmetric Bragg; Entry 85 addendum) are T_flat 2 / T_half 5 / T_90 7 ns, inside one transit; the
beam condition is 0.1–0.4 µrad for both the single-ray factor and the centre law. It is the second experiment, not
the first.

## 5. Decision for Isaac

Regime A alone is a short, low-risk measurement that settles the practical question for
picosecond ultrasonics (no crystal factor with an ordinary beam) and tests Entry 82. Regime B is
the paper's crystal factor and needs an IXS-type beamline; propose it as the second step, 8 keV,
with the ±ω ratio as its primary observable (largest effect, ratio measurement, no absolute
calibration).
