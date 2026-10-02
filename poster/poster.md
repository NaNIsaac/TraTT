---
title: Transient Formation Dynamics of Acoustic Bragg Satellites in Dynamical X-ray Diffraction
authors: **Isaac Yong**¹, Ouyoung Kwon¹, Hyunjung Kim¹*
affiliation: ¹Department of Physics, Sogang University, Seoul, Korea  ·  *corresponding author
footer: P2.17 · 2026 KPS Fall Meeting, Jeju, October 22, 2026 · This work was supported by the National Research Foundation of Korea grant RS-2021-NR059920.
---

<!--
Poster manuscript. This file is the only source of printed text; poster/build.py reads it and
never edits it. Every item carries {E..} tags pointing at notes/lab_notebook.md (checked).
Syntax: see the docstring of poster/build.py. Poster figures: python explore/make_poster_figs.py
-->

@wide

## Overview: how do acoustic satellites form?

A travelling acoustic wave enters a perfect crystal. The steady-state satellites are known, but how do they **form** while the wavefront crosses the crystal? Is the transit time $t_{tr}$ the only clock, or does dynamical diffraction add one? {E0}

![](figures/poster/f0_scheme.png){E0}

> **Answer.** Behind the front, satellite $n$ sees a crystal of the filled depth $v_s t$ with coupling $\pi\alpha|J_n(S)|$. The transit time ends the transient; it does not shape it. Inside the transit the clock is dynamical, $t_{tr}/(\alpha|J_n(S)|)$. {E1}{E9}{E33}

## Model

- Stroboscopic two-beam Takagi–Taupin, Laue and Bragg, with $\phi = S\,A(t-\zeta)\sin(q\zeta-qt)$; $\alpha = L/\Lambda_P$, $S = h\cdot u_0$, $R = \Lambda_{ac}/\Lambda_P$. {E0}
- Checked against analytic limits ($10^{-9}$) and an independent exact Floquet solver. {E0}{E42}

## Main result: the filled-depth mapping

- Laue: $I(nq,t) = \sin^2(\pi\alpha_n t)$, a temporal Pendellösung with period $\Lambda_P/(v_s|J_n(S)|)$, independent of $L$. {E1}
- Bragg: $I(nq,t) = \tanh^2(\pi\alpha_n t)$; the naive $t_f \sim t_{tr}/\alpha$ holds with $\alpha \to \alpha|J_n(S)|$. {E9}

![Points: full Takagi–Taupin for four $(S, \alpha, n)$; lines: closed forms, no fit parameter.](figures/poster/f1_collapse.png){E1}{E9}

## Formation is not peak growth

- At the dark point ($\alpha_n t = 1$) the centre falls to 0.1 % of its maximum, yet ~52 % of the integrated satellite remains in two side lobes: Pendellösung fringes flow inward. {E23}

![Satellite $n=1$, $S=1.5$, $\alpha=3$: the line narrows and blinks as the front fills the crystal.](figures/poster/f2_fringes.png){E22}{E23}

@column

## Formation remembers the path

Satellite $n$ reads the last $t_{tr}$ of the drive, $\Theta_n = \int_0^1 J_n(S\,A(t-\zeta))\,d\zeta$ ($\sin^2$ Laue, $\tanh^2$ Bragg). An order forbidden at the final amplitude, $J_1(S_f)=0$, is the brightest line during formation. {E28}{E30}{E34}

![Laue, $\alpha=0.5$: points full TT, lines the $\Theta_n$ law.](figures/poster/f3_memory.png){E28}{E30}

## Where the independent-satellite picture ends

- Exact phase-averaged steady state: channels stay independent (error below 5%) for $R$ below 0.34–0.57, depending on $S$; all results above sit at $R \lesssim 0.1$. {E53}{E54}
- $R \sim 1$ is the known X-ray–acoustic resonance (Entin; Iolin & Entin). {E50}

![](figures/poster/f4_regime.png){E54}

## Not new features, and toward an experiment

- Peak shifts and $\pm n$ asymmetry are tail interference; acoustic phase walls are stacking faults with phase $n\,\delta\phi$ (Authier & Simon 1968, within 1–3 %). {E14}{E56}
- Si 440 at 17–30 keV: ~0.6–2 GHz drive, $\tau_1$ = 17–30 ns; $S$ calibrates in situ from $\tau_2/\tau_1 = |J_1/J_2|$ to ~1 %. {E13}{E20}
