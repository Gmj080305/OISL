# Report revisions (log)

**Status: all items below were applied in v2.5 (November 2026).** The report now
matches the code. This file records what was wrong and why, so the reasoning is
kept. Line numbers refer to the October 2026 version of the chapters.

Two fixes went further than the original suggestion:
- PAPR is now the 10⁻⁴ CCDF point (7.5 / 7.4 / 6.4 dB), computed in the code.
- The heuristic reconstruction penalty is excluded from all results
  (`P.recon = false`) and described in the report as a failed closed-form attempt.

## A. Statements contradicted by the simulations (fix first)

| File:line | Current text | Problem | Suggested direction |
|---|---|---|---|
| ch3:36 | "the noise model of Section 3.1.3 is accurate for κ ≥ 9 dB" | ζ = 1 + 2/κ verified within 0.15 dB only for κ ≥ 12 dB; 0.8 dB off at 9 dB | State κ ≥ 12 dB |
| ch3:36, 101, Table 3.2 (209) | PAPR 8.8 / 8.7 / 8.5 dB; recon SNR 61.8 dB | PAPR is the peak of one random waveform and depends on the generator and seed (Octave: 7.9 / 7.9 / 6.8 dB). The QPSK recon SNR depends on it (62.0 dB in the current run) | Quote PAPR as seed-dependent, or report a percentile (e.g. 10⁻⁴ CCDF point) instead of the maximum |
| ch3:97–101, ch1:34, ch4:10 | Reconstruction penalty "determines the optimal CSPR"; κ* ≈ 9–10 dB; listed as a main contribution | The waveform simulation shows the formula is 35–60 dB too pessimistic at low κ. It also lowers BPSK to 29.2 dB on-axis, while the waveform gives 31.05 dB (visible in old Fig 3.5b). Waveform optimum: κ* ≈ 7 dB at A₀ and 0.1·A₀, 9 dB at 0.03·A₀ | Present the formula as a heuristic upper bound on the penalty; take κ* from the waveform sweep (`RUN.cspr_wave`); remove it from the contribution list or reframe it. Wen et al. measure 7–8 dB, which supports the waveform value |
| ch3:101 | BPSK "field crosses zero when s < −√κ" | The formula uses \|s\| > √κ, and the SSB shift makes the signal complex | Align the wording with the formula |
| ch3:120, 261 | Bell curve and κ* ≈ 9–10 dB "for the default parameters"; refers to a plot that is not in the report | No figure shows it | Add the waveform CSPR figure (`fig_opt_snr_vs_cspr_waveform`) or drop the claim |
| ch3:257, ch4:20 | Waveform MC "agrees within 0.2 dB between A₀ and 0.01·A₀" | ≤ 0.1 dB from A₀ to 0.03·A₀; 0.4 dB at 0.01·A₀ (QPSK) | Correct the numbers; say the 0.01·A₀ gap is KK's nonlinearity at low SNR |
| ch3:272, 315 | B1: "transmitter sends only the SSB data signal at full power" | Old code sent P_T/(1+κ) and amplified the LO through the EDFA. The new code matches the text: data at full P_T, LO 1.77 mW (CSPR 12 dB at A₀), no EDFA, ideal LO | Keep the text; add LO power and "ideal LO (no offset, no linewidth)" |
| ch1:36, ch3:276, ch4:12 | B2 = "data-aided phase rotation after KK" | The final gain/phase estimate already removes a common rotation, so B2 was numerically identical to A1/B3 (both −36.3 dB). B2 is now "the carrier level is added digitally after detection, before KK" | Rewrite B2; it still fails (BER ≈ 0.5), for the reason the text already gives |
| ch3:267–289 | Fig 3.6 / 3.7 numbers (A2 33.8 dB etc. from the old Python runs) | MATLAB chain: A2 31.1 dB (matches the analytic value), B1 37.3 dB, A1/B3 −27.9 dB, B2 −24.4 dB; σ_φ 0.020 / 0.009 / 1.82 / 1.83 rad; metrics from 7680 payload symbols after a 512-symbol pilot estimate | Update numbers and captions. Point out that B1's +6.2 dB over A2 is the CSPR penalty made visible, with an ideal LO |
| ch3:113 | "With the pre-amplifier, the ASE term (∝ y⁻¹) dominates" | True only above y ≈ 0.02·A₀. Below that the thermal term (∝ y⁻²) dominates even with the EDFA (Fig 3.4), which is the region where outages happen | Qualify: ASE-limited near the axis, thermal-limited in deep fades |
| ch3:373 | "trades 6–10 dB of transmit power" | Not derived anywhere | Use the simulated 6.2 dB (A2 vs B1 at A₀, ideal LO), or derive the number |

## B. Wording and smaller corrections

| File:line | Issue | Fix |
|---|---|---|
| ch1:38, ch3:293, ch4:14 | "94 % of the maximum possible rate" | It is 94 % of the rate at w_d = 0.10 m, the largest value in the sweep; the rate keeps rising with w_d |
| ch3:182, ch3:293, ch4:14 | w_d "≈30 % below beam waist" / "20–50 % below the beam waist w_z" | 0.07 m is 30 % below the **Tx waist w₀** = 0.1 m; w_z is 4.87 m |
| ch1:38, ch3:306, ch4:16 | "less than 2 %" (caption: "across the power range") | 13.18 → 13.45 Gb/s is 2.05 % at 0.3 W; at 0.1 W it is 11.49 → 11.75 Gb/s (2.3 %) | Write "about 2 %" |
| ch1:40 | "MATLAB + Python mirror", "zero external dependencies" | The Python mirror is not in the repository. The Communications Toolbox is optional (one extra plot). Write "MATLAB/Octave, no toolboxes required" |
| ch3:25 | Residual offset ≤ 266 kHz "tracked digitally" | Not simulated. Either model it or say it is assumed (see `docs/open_questions.md`, item 4) |
| ch3:183–191 | B_f = 30 GHz has no justification | It follows from Doppler: B_o + 2 × 9.9 GHz ≈ 30.3 GHz. One sentence |
| ch2:31 | "α ∈ [0.1, 5]" for LEO at 1000 km | No citation; cite or remove |
| ch4 Limitations | Missing items | Add: ζ validity (κ ≥ 12 dB), PAPR seed dependence, BER floor with 2¹³ symbols (~10⁻⁴), Doppler time-scaling not modelled |
| ch4 Future work | Missing items | Add: OptiSystem plan (`docs/open_questions.md`, Next steps); CKK / GKK / MP-KK once real hardware impairments are present (Wen et al. 2023) |
| main.tex bibliography | Wen et al. not listed | Add if cited (entry in `references/README.md`) |

## C. Figures

The figures the LaTeX uses are now in `report/figs/`, copied from `figures/`
under the names below. To refresh them, run the script with `PLOT.save = true`
and copy again.

| LaTeX file | Generated by the current code | Notes |
|---|---|---|
| figs/fig31.pdf | fig31_rate_vs_freq_BPSK.pdf | `RUN.freq` |
| figs/fig32.pdf | fig32_rate_vs_freq_QPSK.pdf | `RUN.freq` |
| figs/fig33.pdf | fig33_rate_vs_snr.pdf | `RUN.snr` |
| figs/fig34.pdf | fig34_snr_budget_vs_y.pdf | `RUN.ybudget` |
| figs/fig35a.pdf, fig35b.pdf | fig35a_channel_mc_cdf.pdf, fig35b_waveform_mc.pdf | two subfigures |
| figs/fig36.pdf | fig36_constellation_study.pdf | one figure (old Fig 3.7 merged into it). Black background |
| figs/fig38.pdf | fig38_rate_vs_wd.pdf | `RUN.wd` |
| figs/fig39.pdf | fig39_rate_vs_PT_cband.pdf | `RUN.itu` |
| figs/fig_cspr_wave.pdf | fig_opt_snr_vs_cspr_waveform.pdf | new: CSPR optimum |
| (new, optional) | fig310–fig313 | only if Section 10 (MP-KK) is included |
