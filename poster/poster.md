---
title: Transient Formation Dynamics of Acoustic Bragg Satellites // in Dynamical X-ray Diffraction
authors: **Isaac Yong**¹, Ouyoung Kwon¹, Hyunjung Kim¹*
affiliation: ¹Center for Ultrafast Phase Transformation, Department of Physics, Sogang University, Seoul, Korea  ·  *corresponding author
footer: P2.17 · 2026 KPS Fall Meeting, Jeju, October 22, 2026 · This work was supported by the National Research Foundation of Korea grant RS-2021-NR059920.
---

<!--
Poster manuscript. This file is the only source of printed text; poster/build.py reads it and
never edits it. Every item carries {E..} tags pointing at notes/lab_notebook.md (checked).
Syntax: see the docstring of poster/build.py. Poster figures: python explore/make_poster_figs.py

Three claims, in this order, and everything else is evidence for them:
  1 each satellite has its own dynamical clock;
  2 the transient remembers the drive history;
  3 the independent-satellite picture holds up to R ~ 0.3.
-->

@wide

## Overview: a new clock for acoustic satellites

A travelling acoustic wave enters a perfect crystal. The steady-state satellites are known; here we resolve how they **form** while the wavefront crosses the crystal, and identify the clock that governs it. {E0}

Temporal Pendellösung of the main reflection: DeCamp et al., Nature 413, 825 (2001). Steady-state acoustic satellites: Entin; Zolotoyabko et al. {E7}

- Two-beam Takagi–Taupin, travelling acoustic deformation, Laue and Bragg; $\alpha = L/\Lambda_P$, $S = h\cdot u_0$, $R = \Lambda_{ac}/\Lambda_P$. The lattice is static while the field is inside, $\beta \leq 1.2\times10^{-2}$. {E0}{E62}
- Verified to $10^{-9}$ against analytic limits and $3\times10^{-7}$ against an independent exact Floquet solver. {E0}{E42}

![Snapshot at $\tilde t$ = 0.62: the profile inside the slab is the drive, the curve on the right the solved $I(\eta)$ of that same instant ($\alpha$ = 0.6, $S$ = 1.5, $\alpha_{ac}$ = 5).](figures/poster/f0_scheme.png){E0}

> **Key result.** Each satellite diffracts from the filled depth $v_s t$, on its own clock $\tau_n = \Lambda_P/(v_s\left|J_n(S)\right|)$. **Transit time sets the window; extinction length sets the clock.** {E1}{E9}{E33}

## Filled-depth mapping: one law, no fit

- Behind the front, satellite $n$ is a crystal of filled depth $v_s t$ with coupling $\alpha_n = \alpha\left|J_n(S)\right|$: Laue reads $I_n(t) = \sin^2(\pi\alpha_n t)$, Bragg reads $\tanh^2(\pi\alpha_n t)$. {E1}{E9}

![Points: full Takagi–Taupin for four $(S, \alpha, n)$; lines: the closed forms, **no fit parameter**. The Laue period $\Lambda_P/(v_s\left|J_n(S)\right|)$ does not depend on $L$. Formation order follows the Bessel ladder, so $S$ reorders it — at $S=3$: $n=2$, then 1, then 3.](figures/poster/f1_collapse.png){E1}{E9}{E24}

## Formation is not peak growth

- At the dark point ($\alpha_n t = 1$) the centre falls to **0.1–0.2 %** of maximum, yet **≈52 %** of the integrated satellite sits in two side lobes: fringes flow inward. {E23}

![Satellite $n=1$, $S=1.5$, $\alpha=3$: the line narrows and blinks as the front fills the crystal.](figures/poster/f2_fringes.png){E22}{E23}

## When does the independent-satellite picture fail?

- Accurate to 5 % up to $R \lesssim 0.3$; computed boundary **0.34–0.57**, depending on $S$. {E53}{E54}
- Past it the channels talk through sidebands of the forward beam — the same third-order path, $\pi\alpha (R/2)^2 C_n(S)$, that brings a Bessel-dark order back. {E44}{E45}{E50}{E52}

![Population error of the independent-channel picture. Each point is one $(\alpha, \alpha_{ac}, S)$, exact acoustic-phase average.](figures/poster/f4_regime.png){E54}

@column

## Formation remembers the path

> **Same final amplitude $\neq$ same diffraction state.** {E30}

- Satellite $n$ reads the last $t_{tr}$ of the drive, $\Theta_n = \int J_n(S\,A(t-\zeta))\,d\zeta$: an order forbidden at the final amplitude ($J_1(S_f)=0$) is the **brightest** line during formation. {E28}{E30}{E34}

![Laue, $\alpha=0.5$: points full TT, lines the $\Theta_n$ law.](figures/poster/f3_memory.png){E28}{E30}

## Outlook: Bragg oscillations as a phonon spectrometer

- The Fourier weight of $R(t)$ is the strain spectrum times a crystal transfer function $H(\omega)$, closed form at first order in $S$; full TT agrees to 2–5 %. {E59}
- It holds inside the Darwin width, where the plane-wave perturbation of Sondhauss & Wark (2003) is off by ×3–8. Si 400 at 17 keV: **×2** bias at 20 GHz. {E64}{E65}

![(a) Fourier weight over the kinematic law vs phonon frequency, $\kappa/2\pi\alpha$ = 0.12: three packets, one curve; divide it out to recover the spectrum. (b) Satellite weight of the first-order forms over full TT, single mode.](figures/poster/f5_outlook.png){E59}{E64}

## Take home

%% $\tau_n \propto 1/\left\vert J_n\right\vert$ | clock: set by each order's own extinction length {E1}{E9}
%% $\Theta_n$ | memory: $\Theta_n=\int J_n(S A)\,d\zeta$ reads the last $t_{tr}$ of the drive {E28}{E30}
%% $R \lesssim 0.3$ | validity: 5 % accuracy for independent satellites {E54}

A travelling acoustic front makes every Bragg satellite a dynamical probe: its own clock, its own memory of the drive, valid while the channels stay resolved. {E33}

- Si 440 at 17–30 keV, 1–2 GHz: $\tau_1$ = 17–30 ns. {E20}
