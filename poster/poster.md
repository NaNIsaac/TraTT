---
title: Transient Formation Dynamics of Acoustic Bragg Satellites in Dynamical X-ray Diffraction
authors: AUTHORS (to fill)
affiliation: AFFILIATION (to fill)
contact: CONTACT (to fill)
---

<!--
Poster manuscript. This file is the only source of printed text; poster/build.py reads it and
never edits it. Rules (checked by the build):
- every bullet and paragraph carries at least one {E..} tag pointing at notes/lab_notebook.md;
  tags are stripped from the print and listed in the build log;
- figures: ![caption](figures/file.png){E..}; the file must exist;
- "## n | Title" starts a panel; "@column" forces the next panel into the next column.
Math: $...$ (matplotlib mathtext).
-->

## 1 | Question

A travelling acoustic wave enters a perfect crystal. The steady-state satellites are known. How do they **form** while the wavefront crosses the crystal: is the transit time $t_{tr}$ the only clock, or does dynamical diffraction add one? {E0}

## 2 | Model

- Stroboscopic two-beam Takagi–Taupin with the phase grating $\phi = S\,A(t-\zeta)\sin(q\zeta - qt)$; Laue and Bragg. {E0}
- Dimensionless: $\alpha = L/\Lambda_P$, $S = h\cdot u_0$, $\alpha_{ac} = L/\Lambda_{ac}$, $R = \Lambda_{ac}/\Lambda_P$. {E0}
- Verified against analytic limits to $10^{-9}$, 4th-order convergence, and an independent exact Floquet solver for the filled crystal. {E0}{E42}

## 3 | Main result: the filled-depth mapping

Behind the front, satellite $n$ sees a crystal of thickness $v_s t$ with coupling $\pi\alpha|J_n(S)|$, and replays its static thickness curve in time. {E1}{E9}

- Laue: $I(nq,t) = \sin^2(\pi\alpha_n t)$, a temporal Pendellösung with period $\tau_n = \Lambda_P/(v_s|J_n(S)|)$, independent of $L$. {E1}
- Bragg: $I(nq,t) = \tanh^2(\pi\alpha_n t)$; the naive $t_f \sim t_{tr}/\alpha$ holds with $\alpha \to \alpha|J_n(S)|$. {E9}
- No fit parameter: all $(\alpha, S, n)$ collapse onto one curve. {E1}{E9}

![Satellite intensity at the channel centre vs time; full TT and the filled-depth law.](figures/poster_mapping.png){E1}{E9}

## 4 | Formation is not peak growth

- At the dark point ($\alpha_n t = 1$) the centre is at 0.1 % of its maximum, yet ~52 % of the integrated satellite remains, in two side lobes: Pendellösung fringes flow inward. {E23}
- One front, several stages: at $S = 3$ the $n=2$ line is split around a dark centre while $n=3$ peaks and $n=1$ decays. {E24}

![Time-resolved rocking curve $I(\eta,t)$, full TT next to the closed three-mode model.](figures/poster_map.png){E22}{E23}

@column

## 5 | Formation remembers the path

The crystal is a record of the last $t_{tr}$ of the drive. Satellite $n$ reads $\Theta_n = \int_0^1 J_n(S\,A(t-\zeta))\,d\zeta$: $\sin^2\Theta$ in Laue, $\tanh^2\Theta$ in Bragg, exact at resonance. {E28}{E30}

- The steady-state pattern does not determine the transient one: an order forbidden at the final amplitude ($J_1(S_f) = 0$) is the brightest line during formation. {E28}{E34}
- Same strain, different pattern: $I_1 = 0.010$ on the way up vs 0.494 on the way down a triangular pulse. {E30}
- A $\pi$-coded packet cancels $n=1$ while both halves are inside; $n=2$ peaks, then $n=1$ revives. {E30}

![Transient-only satellite, up/down asymmetry and the $\pi$-coded packet.](figures/poster_memory.png){E28}{E30}

## 6 | Where the independent-satellite picture ends

Steady state, filled crystal, exact acoustic-phase average $\Sigma_m |b_m|^2$. {E43}{E53}{E54}

- Independent channels hold (population error below 5%) for $R$ below 0.34–0.57, depending on $S$. All formation results above sit at $R \lesssim 0.1$. {E54}
- $R \sim 1$ is the known X-ray/neutron–acoustic resonance (Entin; Iolin & Entin). {E50}
- At a Bessel zero, the dark order reappears through a third-order indirect coupling $g_{eff} = \pi\alpha (R/2)^2 C_n(S)$, which also shifts the zero. {E45}{E52}

![Population error of the independent-channel picture vs $R$ (exact phase average where computed).](figures/regime_phase_diagram.png){E54}

@column

## 7 | What looked like features but are not

- Peak shifts, $\pm n$ asymmetry and the acoustic-period ripple are one kinematic effect: interference with the main line's off-Bragg tail. {E14}
- The single-phase "breakdown" of order additivity is micromotion; it averages away. {E43}{E53}
- Chirps and acoustic phase walls only relieve extinction by dephasing; a phase wall is a stacking fault with phase $n\,\delta\phi$ (Authier & Simon 1968 reproduces our TT to 1–3 %). {E41}{E56}

## 8 | Toward an experiment

- Resolved channels need $\Lambda_{ac} \lesssim 0.05$–$0.1\,\Lambda_P$: the drive frequency is set by the Pendellösung length, not the thickness. {E19}{E26}
- Si 440 at 17–30 keV: ~0.6–2 GHz, strain ~$1.5$–$3\times10^{-5}$, $\tau_1$ = 17–30 ns, absorption negligible. {E20}
- $S$ calibrates in situ from $\tau_2/\tau_1 = |J_1/J_2|$ to ~1 %, without detector gain. {E13}

## 9 | Answer

The transit time ends the transient; it does not shape it. Inside the transit the clock is dynamical, $t_{tr}/(\alpha|J_n(S)|)$: an oscillation in Laue, a saturation in Bragg. The static extinction length of each satellite becomes a time. {E33}

## 10 | Outlook

When does the Fourier weight of a time-resolved reflectivity equal the strain spectrum? Comparison with Chalise et al., arXiv:2603.28071 (2026), in progress. {E57}
