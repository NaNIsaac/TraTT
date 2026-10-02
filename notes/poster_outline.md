# KPS 2026 Fall poster — outline draft (2026-09-30, updated 2026-10-02 for Entries 34–57)

Title (abstract): *Transient Formation Dynamics of Acoustic Bragg Satellites in
Dynamical X-ray Diffraction.* Every claim below points to a notebook entry
(`notes/lab_notebook.md`); nothing here goes on the poster without its entry.

## Panel 1 — Question
A travelling acoustic wave enters a perfect crystal. The steady-state
satellites are known. How do they *form* while the wavefront crosses the
crystal — is the transit time t_tr the only clock, or does dynamical
diffraction add one?

## Panel 2 — Model
Stroboscopic two-beam Takagi–Taupin; phase grating
phi = S A(t - zeta) sin(q zeta - q t); dimensionless alpha = L / Lambda_P,
S = h.u0, alpha_ac = L / Lambda_ac; Laue and Bragg; analytic ground truth to
1e-9, 4th-order convergence (C-1.2). Resolved channels, alpha / alpha_ac <= 0.1.
Independent check for the filled uniform wave: an exact Floquet solver
(`src/tratt/floquet.py`, Entry 42). The acoustic-phase average is exactly
sum_m |b_m|^2 (Entry 53, test).

## Panel 3 — Main result: the filled-depth mapping (figure: `figures/poster_mapping.png` a, b)
Behind the front, satellite n sees a crystal of thickness v_s t with coupling
pi alpha |J_n(S)|. It replays its static thickness curve in time:

- Laue: I(nq, t) = sin^2(pi alpha_n t) — temporal Pendellösung, period
  tau_n = Lambda_P / (v_s |J_n(S)|), independent of L. (C-1.10, Entry 1)
- Bragg: I(nq, t) = tanh^2(pi alpha_n t); formation time 0.95 t_tr / alpha_n.
  The naive t_f ~ t_tr / alpha is right with alpha -> alpha |J_n(S)|. (Entry 9)

No fit parameter; all (alpha, S, n) collapse onto one curve.

## Panel 4 — Consequences
- **Formation is not peak growth.** At the dark point (alpha_n t = 1) the
  centre is at 0.1 % of its maximum while ~52 % of the integrated satellite
  remains, in two side lobes. The line shape is a conveyor of Pendellösung
  fringes flowing inward; each centre maximum is a fringe arrival. (Entry 23;
  `figures/poster_map.png`)
- **One front, several stages.** At S = 3 the n = 2 line is split around a dark
  centre while n = 3 peaks and n = 1 decays — tau_n made visible.
  (Entry 24; `figures/poster_ladder_map.png`, replaces the formation-order plot)
- Main line: rotations add, sin^2(pi alpha [1 - (1 - J_0) t]) — CW analogue of
  the DeCamp et al. (2001) dividing layer. (Entry 1)
- Width: Laue narrows to 0.71 of the finite-depth limit; Bragg stops at the
  satellite Darwin width ~13.4 alpha_n. (Entries 4, 17)
- Finite rise: sin^2(pi alpha |int J_n(S A) dzeta|); a ramp through the Bessel
  maximum gives a 1.5x overshoot at weak coupling, through a zero a dark
  window. (Entries 5, 15)
- At a Bessel zero the slow clock stops; only the acoustic-period ripple is
  left. (Entry 25)
- Regime map of the formation law: accuracy ~ 0.9 Lambda_ac / Lambda_P, independent
  of the clock. (Entry 26; `figures/regime_map.png`)

## Panel 4b — Formation remembers the path (figure `figures/poster_memory.png`)
The crystal is a record of the last t_tr of the drive; satellite n reads
Theta_n = int_0^1 J_n(S A(t - zeta)) dzeta, sin^2 in Laue, tanh^2 in Bragg (exact).
- Transient-only satellite: S_f at a J_1 zero — absent before and after, bright
  in between; and it is the brightest order during the transient. (Entry 28)
- Same strain, different pattern: I_1 = 0.010 on the way up vs 0.494 on the way
  down a triangular pulse. (Entry 30a)
- Pi-coded packet: n = 1 cancels while both halves are inside, n = 2 peaks,
  then n = 1 revives. (Entry 30d)
- Brightest order switches 1 -> 2 -> 3 -> 4 -> 5 during a ramp, t_tr / 2 late.
  (Entry 31)
All at alpha = 0.5 (no Pendellösung) and identical in Bragg.
Framing (Entry 34): claim "the steady-state pattern does not determine the
transient one". Do not claim a new kind of peak; the slow-ramp limit is the
adiabatic Bessel ladder. Robust within the model: a >30:1 transient/final
contrast survives a 5 % error in S or 10 % acoustic loss (Entries 34, 35).
Two-domain cancellation (Entry 37) darkens only the centre; the eta-integral
obeys Parseval, so integrated darkness needs a local J_n = 0.

## Panel 4c — Where the independent-satellite picture ends (figure `figures/regime_phase_diagram.png`)
Steady state, filled crystal, exact phase average (Entries 42–54).
- Order-window population error E_pop < 0.05 (independent channels hold) for
  R = Lambda_ac / Lambda_P below 0.44 (S = 0.5), 0.34 (S = 1), 0.44 (S = 2),
  0.34 (S = 3), 0.57 (S = 5, 8); bins x1.3, S-dependent, not a universal law
  (Entry 54 table). All formation results above sit at R <~ 0.05–0.1.
- R ~ 1 is the known X-ray/neutron–acoustic resonance (Entin; Iolin–Entin;
  Raitman et al.; Entry 50). Cite, do not claim.
- First failure at a Bessel zero: the dark order reappears through a
  third-order indirect coupling g_eff = pi alpha (R/2)^2 C_n(S), with a
  renormalised zero S* - j = (R/2)^2 C_n / J'_{-n} (Entries 45, 52; "possibly new"
  only after reading Entin / Iolin–Entin in full).
  Observability limits (Entries 53, 54; alpha = 1, S = j_11):
  - it dominates the phase-averaged exit intensity only within
    |alpha_ac - integer| <~ 1.3e-3 (90 %) and an angular half-width
    w/q <~ 2e-4;
  - in the |eta - q| < 0.45 q integral its share is 4.4e-5.
  Present it as theory, not as an experiment proposal.

## Panel 5 — What looked like features but are not
Peak shifts, ±n asymmetry, acoustic-period ripple and the +12 % integral
excess are one kinematic mechanism: interference with the main line's
off-Bragg tail, reproduced by a closed three-mode model and vanishing as
alpha / alpha_ac -> 0. (Entries 2, 3, 8, 14)
- The single-acoustic-phase "breakdown" of order additivity is micromotion; it
  averages away. The old 4-phase average was not converged (aliasing of
  m - m' = 0 mod 4) and overstated breakdown at small alpha_ac and large S
  (Entries 43, 53).
- Chirps, q-walls and phase walls (Entries 38–41) are one known mechanism:
  dephasing relieves extinction. The phase wall is quantitatively a stacking
  fault with fault phase n dphi: Authier & Simon (1968) eq. (4-5) reproduces
  the TT integrals to 1–3 % (Entries 55, 56). Order-space localisation is
  absent (Entry 51).

## Panel 6 — Experiment
Constants from X0h (Entry 20). The resolved-channel rule is
Lambda_ac <= 0.05–0.1 Lambda_P, so the drive frequency is set by the
Pendellösung length, not the thickness. Si 220 at 8 keV needs 4–8 GHz; **Si 440
at 17–30 keV needs ~0.6–2 GHz**, strain ~1.5–3e-5, tau_1 = 17–30 ns, absorption
negligible — the recommended example. S calibrated in situ from
tau_2 / tau_1 = |J_1 / J_2| to ~1 %, gain-free (Entry 13).
Figure: `figures/poster_map.png` — I(eta, t) from full TT next to the
closed three-mode model (Entry 22).

## Panel 7 — Answer to the abstract's question
(Consider the story of notebook Entry 33 as the lead message.)
The transit time ends the transient; it does not shape it. Inside the transit
the clock is dynamical: t_tr / (alpha |J_n(S)|) — oscillation in Laue,
saturation in Bragg. The static extinction length of each satellite, known
from the steady state, becomes a time.

## Outlook (one line on the poster)
Reading acoustic spectra from time-resolved diffraction: when does the Fourier
weight of R(t) equal the strain spectrum? Compare with Chalise et al.,
arXiv:2603.28071 (2026) (Entries 57–58).

## Before printing
- Read DeCamp 2001 and the ultrasonic-superlattice literature in full; fix the
  novelty sentence (Entry 7).
- X-ray constants verified with X0h (Entry 20); c44 and GHz shear-wave launch through 0.2–0.3 mm Si still to check.
- Read Polikarpov, phys. stat. sol. (b) 143 (1987) / 149 (1988) (paywalled; needs library access) — Entry 21.
- Decide whether Bragg (deprioritised) gets a full panel or a sub-panel.
- Read Entin (JETP Lett. 26, 392 (1977); JETP 77, 214 (1979)) and Iolin–Entin
  (JETP 85, 1692 (1983)) before any "possibly new" wording in Panel 4c.
- Read Authier & Simon part II (moiré, different h) for the q-wall verdict.
- Panel 4c figure: `figures/regime_phase_diagram.png` now uses the best phase
  average per point (exact > 16 > 4 phases).
