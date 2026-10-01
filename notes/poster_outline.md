# KPS 2026 Fall poster — outline draft (2026-09-30)

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

## Panel 3 — Main result: the filled-depth mapping (figure: `figures/poster_mapping.png` a, b)
Behind the front, satellite n sees a crystal of thickness v_s t with coupling
pi alpha |J_n(S)|. It replays its static thickness curve in time:

- Laue: I(nq, t) = sin^2(pi alpha_n t) — temporal Pendellösung, period
  tau_n = Lambda_P / (v_s |J_n(S)|), independent of L. (C-1.10, Entry 1)
- Bragg: I(nq, t) = tanh^2(pi alpha_n t); formation time 0.95 t_tr / alpha_n.
  The naive t_f ~ t_tr / alpha is right with alpha -> alpha |J_n(S)|. (Entry 9)

No fit parameter; all (alpha, S, n) collapse onto one curve.

## Panel 4 — Consequences
- Formation order follows the Bessel ladder, reorderable by S
  (S = 3: n = 2, 1, 3). (Entry 1; figure c)
- Main line: rotations add, sin^2(pi alpha [1 - (1 - J_0) t]) — CW analogue of
  the DeCamp et al. (2001) dividing layer. (Entry 1)
- Width: Laue narrows to 0.71 of the finite-depth limit at alpha_n t ~ 0.86;
  Bragg stops narrowing at the satellite Darwin width ~13.4 alpha_n. (Entries 4, 17)
- Finite rise: sin^2(pi alpha |int J_n(S A) dzeta|) for all channels; a ramp
  through the Bessel maximum gives a 1.5x overshoot even at weak coupling,
  through the zero a dark window. (Entries 5, 15; figure d)

## Panel 5 — What looked like features but are not
Peak shifts, ±n asymmetry, acoustic-period ripple and the +12 % integral
excess are one kinematic mechanism: interference with the main line's
off-Bragg tail, reproduced by a closed three-mode model and vanishing as
alpha / alpha_ac -> 0. (Entries 2, 3, 8, 14)

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
The transit time ends the transient; it does not shape it. Inside the transit
the clock is dynamical: t_tr / (alpha |J_n(S)|) — oscillation in Laue,
saturation in Bragg. The static extinction length of each satellite, known
from the steady state, becomes a time.

## Before printing
- Read DeCamp 2001 and the ultrasonic-superlattice literature in full; fix the
  novelty sentence (Entry 7).
- X-ray constants verified with X0h (Entry 20); c44 and GHz shear-wave launch through 0.2–0.3 mm Si still to check.
- Read Polikarpov, phys. stat. sol. (b) 143 (1987) / 149 (1988) (paywalled; needs library access) — Entry 21.
- Decide whether Bragg (deprioritised) gets a full panel or a sub-panel.
