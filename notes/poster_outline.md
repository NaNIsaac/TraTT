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
Si 220, L = 50–200 um, shear wave ~2 GHz, strain ~1e-4: tau_1 ~ 5 ns,
storage-ring timing suffices; Borrmann absorption keeps contrast (Entry 11b).
S calibrated in situ from tau_2 / tau_1 = |J_1 / J_2| to ~1 %, gain-free
(Entry 13). Material constants still to be verified (X0h/XOP).

## Panel 7 — Answer to the abstract's question
The transit time ends the transient; it does not shape it. Inside the transit
the clock is dynamical: t_tr / (alpha |J_n(S)|) — oscillation in Laue,
saturation in Bragg. The static extinction length of each satellite, known
from the steady state, becomes a time.

## Before printing
- Read DeCamp 2001 and the ultrasonic-superlattice literature in full; fix the
  novelty sentence (Entry 7).
- Verify Lambda_P, kappa, v_T for the chosen reflection.
- Decide whether Bragg (deprioritised) gets a full panel or a sub-panel.
