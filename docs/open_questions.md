# Open questions, limitations and next steps

## Questions for the supervisor

1. **Role of OptiSystem.** Is the OptiSystem co-simulation the primary result, or
   a validation of the analytic model? This decides which numbers go in the
   abstract.
2. **Operating CSPR.** Keep κ = 12 dB, or move to the waveform-optimal value
   (about 7 dB on-axis, 9 dB at 0.03·A₀, +2–3 dB SNR)? Should the Gaussian
   reconstruction-penalty formula be replaced by a waveform-calibrated table, or
   dropped and the analysis restricted to κ ≥ 10 dB?
3. **Receiver-LO comparison (B1).** Which LO power, and should the LO be ideal or
   have a realistic linewidth and residual offset? Those LO costs are what
   numerically justify choosing the transmitter carrier.
4. **Doppler time-scaling.** At 15 km/s the symbol clock shifts by 50 ppm and a
   residual offset of about 260 kHz remains after KK. Model it? It follows from
   numbers already in the report and is not an arbitrary penalty.
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

- The reconstruction penalty $\mathrm{PAPR}_s\,P_{viol}(\kappa)$ is heuristic. At
  κ = 6 dB the waveform gives 44.8 dB reconstruction-limited SNR where the formula
  gives 10.5 dB. For BPSK it lowers the analytic on-axis SNR to 29.2 dB while the
  waveform gives 31.05 dB.
- $\zeta = 1 + 2/\kappa$ holds within 0.15 dB for κ ≥ 12 dB; it is 0.8 dB off at
  9 dB and 4.7 dB off at 6 dB.
- The measured PAPR is the peak of one random waveform. It depends on the random
  generator and seed (QPSK: 8.5 dB in the old Python runs, 6.8 dB in Octave).
- In the deepest fades (≈0.01·A₀, SNR ≈ 7 dB) KK falls 0.4 dB short of the linear
  noise model.
- With 2¹³ symbols per run, BER below about 10⁻⁴ cannot be measured; lower values
  are SNR-based.
- Not modelled: Doppler time-scaling and residual offset, EDFA gain saturation,
  laser linewidth for the B1 LO (unless set), any hardware bandwidth or
  nonlinearity, atmospheric effects (vacuum link).
- The analytic function `oisl_kk_channel.m` models the transmitter-carrier link
  only.

## Planned next steps

1. Run the code once in MATLAB and compare with `results.md`.
2. Update the report text (`report/REVISIONS.md`).
3. OptiSystem co-simulation: CW laser (linewidth, RIN) → QPSK + RRC → IQ modulator
   with the carrier set by bias → MATLAB component scaling the field by
   $\sqrt{h_{PL}y}$ → EDFA → 30 GHz optical filter → PIN photodiode → Bessel
   low-pass → MATLAB component with KK and DSP. First reproduce 31.1 dB at A₀ with
   every component ideal, then switch impairments on one at a time.
4. Sweep y/A₀ with Parameter Sweep (about 10 points) and weight the results with
   $f_Y$ to get ergodic rates (semi-analytic method).
5. Only with real component impairments present: compare CKK, GKK and MP-KK
   (Section 10 code).
