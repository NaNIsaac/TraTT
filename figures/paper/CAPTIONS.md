# Captions for figures/paper (data file, script, commit per figure)

## fig1_estimator.png
Fig. 1. (a) Reflectivity change of a semi-infinite Si 400 crystal at 17.48 keV after the Entry 66 strain packet (87 GHz, 0.9 cycles, peak strain 1e-5) enters at t = 0, at a detuning of 40 GHz: the signal oscillates at the detuning frequency and its amplitude decays on the absorption/extinction scale (Entry 75; explore/nb/e75_window_cf.py integrate_packet, continuation past the packet at first order for the picture). (b) Fourier weight of R(t) at omega = |eta| divided by the kinematic law, against omega / omega_D, for three packet spectra (R_D = 0.05, 0.2, 0.5) at kappa / 2 pi alpha = 0.12: full TT (points, explore/nb/e58_abs.json, e58_darwin.json) and the first-order closed form (lines, explore/nb/e59_transfer.py predict); Entry 59. Drawn by explore/nb/paper_figs.py fig1 at 63bcb14.

## fig2_planewave_vs_bloch.png
Fig. 2. First-order predictions of the Fourier amplitude of R(t) at the mode frequency W, divided by the full TT value, against the detuning in Darwin half-widths (alpha = 8, kappa = 6, S = 0.02, single mode): the plane-wave form of Sondhauss and Wark (eq. 29), its incident-wave normalised form, and the Bloch-wave (Riccati) form of Entry 59. Left: mode at 5 half-widths (R_D = 0.2); right: 20 half-widths (R_D = 0.05). Only points with |F_TT| above 2 % of its maximum are shown. Data explore/nb/e64_sw.json from explore/nb/e64_sw_compare.py; Entries 64, 69 (S&W transcription, PDF not seen). Drawn by explore/nb/paper_figs.py fig2 at 63bcb14.

## fig4b_formation_collapse.png
Fig. 4b. Formation of the n-th acoustic satellite behind a step front: centre intensity against alpha |J_n(S)| t/t_tr for four (S, alpha, n) cases in Laue (left) and Bragg (right) geometry, full TT (markers) and the independent-satellite laws sin^2 and tanh^2 (line), no fitted parameter. Data explore/nb/poster_f1_points.npz (the cached TT points of the poster figure f1, explore/make_poster_figs.py; Entries 1, 9); the poster figure itself is unchanged. Drawn by explore/nb/paper_figs.py fig4b at 63bcb14.

