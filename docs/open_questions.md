# Open questions, limitations and next steps

## Questions for the supervisor

1. **Role of OptiSystem.** Is the OptiSystem co-simulation the primary result, or
   a validation of the analytic model? This decides which numbers go in the
   abstract.
2. **Operating CSPR.** Keep κ = 12 dB, or move to the waveform-optimal value
   (about 7 dB on-axis, 9 dB at 0.03·A₀, +2–3 dB SNR)? The heuristic
   reconstruction formula has been dropped from the results (v2.5); a
   waveform-calibrated g(y, κ) table would allow analytic rates below 12 dB.
3. **Receiver-LO comparison (B1).** Which LO power, and should the LO be ideal or
   have a realistic linewidth and residual offset? Those LO costs are what
   numerically justify choosing the transmitter carrier.
4. **Jitter model.** Shang et al. model σ_s as a length growing as e^{k₀z/d₀}.
   Real pointing jitter is angular (σ_s = z·σ_θ), which makes α independent of
   range and removes the lower edge of the feasible range window
   (`results.md` §8). Which model should the paper use, and with which σ_θ?
   (Doppler, the fair LO comparison and BER near the FEC threshold are done; see
   `results.md` §6, 7, 9.)
5. **DSP scope.** On this link, GKK, MP-KK, DFE and MLSE give no gain at the
   operating CSPR (no nonlinearity, memory or ISI in the model). Should they be
   evaluated only at the OptiSystem stage with real component models? Is FEC
   (convolutional code with Viterbi decoding) wanted?
6. **Hardware values to fix:** PD/TIA bandwidth, ADC rate and ENOB, EDFA output
   power, laser RIN and linewidth. Is there a reference terminal datasheet?
7. **Jitter model.** Is σ_s0 = 2 m realistic? Pointing jitter dominates the rate
   and the 41 % outage figure.
8. **Metrics.** HD-FEC 3.8×10⁻³ or SD-FEC ~2×10⁻²; EVM, Q-factor? How many
   symbols per point?

## Known limitations of the current model

- No validated closed-form model of KK reconstruction exists in this work. The
  heuristic $\mathrm{PAPR}_s\,P_{viol}(\kappa)$ (off by default since v2.5) gives
  10.5 dB at κ = 6 dB where the waveform gives 44.8 dB.
- $\zeta = 1 + 2/\kappa$ holds within 0.15 dB for κ ≥ 12 dB; it is 0.8 dB off at
  9 dB and 4.7 dB off at 6 dB.
- In the deepest fades (≈0.01·A₀, SNR ≈ 7 dB) KK falls 0.4 dB short of the linear
  noise model.
- With 2¹³ symbols per run, BER below about 10⁻⁴ cannot be measured; the
  dedicated BER runs (§9 of `results.md`) use about 10⁶ symbols.
- Not modelled: Tx/Rx oscillator mismatch (timing recovery), Doppler rate within
  a frame, EDFA gain saturation, hardware bandwidth or nonlinearity, atmospheric
  effects (vacuum link).
- Doppler recovery fails occasionally below the FEC threshold (0.01·A₀,
  SNR ≈ 6.6 dB), where the link is in outage anyway.
- The B1 phase-tracking window is chosen per case with knowledge of the payload
  SNR (optimistic for B1).
- The analytic function `oisl_kk_channel.m` models the transmitter-carrier link
  only.

## Planned next steps

1. Run the code once in MATLAB and compare with `results.md`.
2. Decide the jitter model (item 4) and redo the pointing study with it.
3. A KK vs coherent vs IM/DD comparison under the same fading, with consistent
   noise models for each receiver.
3. OptiSystem co-simulation: CW laser (linewidth, RIN) → QPSK + RRC → IQ modulator
   with the carrier set by bias → MATLAB component scaling the field by
   $\sqrt{h_{PL}y}$ → EDFA → 30 GHz optical filter → PIN photodiode → Bessel
   low-pass → MATLAB component with KK and DSP. First reproduce 31.1 dB at A₀ with
   every component ideal, then switch impairments on one at a time.
4. Sweep y/A₀ with Parameter Sweep (about 10 points) and weight the results with
   $f_Y$ to get ergodic rates (semi-analytic method).
5. Only with real component impairments present: compare CKK, GKK and MP-KK
   (Section 10 code).
