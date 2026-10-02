---
title: Transient Formation Dynamics of Acoustic Bragg Satellites in Dynamical X-ray Diffraction
authors: **Isaac Yong**¹, Ouyoung Kwon¹, Hyunjung Kim¹*
affiliation: ¹Center for Ultrafast Phase Transformation, Department of Physics, Sogang University, Seoul, Korea  ·  *corresponding author
footer: P2.17 · 2026 KPS Fall Meeting, Jeju, October 22, 2026 · This work was supported by the National Research Foundation of Korea grant RS-2021-NR059920.
---

<!--
Poster manuscript. This file is the only source of printed text; poster/build.py reads it and
never edits it. Every item carries {E..} tags pointing at notes/lab_notebook.md (checked).
Syntax: see the docstring of poster/build.py. Poster figures: python explore/make_poster_figs.py
-->

@wide

## Overview: a new clock for acoustic satellites

A travelling acoustic wave enters a perfect crystal. The steady-state satellites are known; here we resolve how they **form** while the wavefront crosses the crystal, and identify the clock that governs it. {E0}

![](figures/poster/f0_scheme.png){E0}

> **Key result.** Behind the front, satellite $n$ is a crystal of the filled depth $v_s t$ with coupling $\pi\alpha|J_n(S)|$. Formation runs on a dynamical clock, $t_{tr}/(\alpha|J_n(S)|)$: the transit time only ends the transient. {E1}{E9}{E33}

## Model

- Time-dependent two-beam Takagi–Taupin, Laue and Bragg, with $\phi = S\,A(t-\zeta)\sin(q\zeta-qt)$; $\alpha = L/\Lambda_P$, $S = h\cdot u_0$, $R = \Lambda_{ac}/\Lambda_P$. {E0}
- Verified to $10^{-9}$ against analytic limits and an independent exact Floquet solver. {E0}{E42}

## Filled-depth mapping: one law, no fit

- Laue: $I(nq,t) = \sin^2(\pi\alpha_n t)$, a temporal Pendellösung with period $\Lambda_P/(v_s|J_n(S)|)$, independent of $L$. {E1}
- Bragg: $I(nq,t) = \tanh^2(\pi\alpha_n t)$, formation time $0.95\,t_{tr}/\alpha_n$. {E9}
- Formation order follows the Bessel ladder: channel $n$ first peaks at $t^*/t_{tr} = 1/(2\alpha|J_n(S)|)$, so $S$ reorders the channels (at $S=3$: $n=2$, then 1, then 3). {E1}{E24}

![Points: full Takagi–Taupin for four $(S, \alpha, n)$; lines: closed forms, no fit parameter.](figures/poster/f1_collapse.png){E1}{E9}

## Formation is not peak growth

- At the dark point ($\alpha_n t = 1$) the centre falls to 0.1 % of its maximum, yet ~52 % of the integrated satellite remains in two side lobes: Pendellösung fringes flow inward. {E23}

![Satellite $n=1$, $S=1.5$, $\alpha=3$: the line narrows and blinks as the front fills the crystal.](figures/poster/f2_fringes.png){E22}{E23}

@column

## Formation remembers the path

Satellite $n$ reads the last $t_{tr}$ of the drive, $\Theta_n = \int J_n(S\,A(t-\zeta))\,d\zeta$. An order forbidden at the final amplitude, $J_1(S_f)=0$, is the brightest line during formation. {E28}{E30}{E34}

![Laue, $\alpha=0.5$: points full TT, lines the $\Theta_n$ law.](figures/poster/f3_memory.png){E28}{E30}

## Validity map of the satellite picture

- With the exact acoustic-phase average, satellites stay independent channels (error below 5%) up to $R$ = 0.34–0.57, depending on $S$. {E53}{E54}
- Near the X-ray–acoustic resonance $R \sim 1$, Bessel-dark orders reappear via a third-order coupling $\pi\alpha (R/2)^2 C_n(S)$. {E45}{E50}{E52}

![](figures/poster/f4_regime.png){E54}

## Conclusion and outlook

- Satellite formation is governed by a dynamical clock set by each order's extinction length, with closed-form laws in Laue and Bragg. {E1}{E9}{E33}
- Acoustic phase steps act as tunable, moving stacking faults (Authier–Simon theory to 1–3 %). {E56}
- Si 440, 17–30 keV, 0.6–2 GHz: $\tau_1$ = 17–30 ns, with in-situ $S$ from $\tau_2/\tau_1$ to ~1 %. {E13}{E20}
- Outlook: a closed-form transfer function turns Bragg oscillations into phonon spectra. {E59}
