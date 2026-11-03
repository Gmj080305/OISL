# Reference results

Numbers produced by the code in `matlab/` (v2.6) with the default parameters
(1530 nm, P_T = 0.3 W, z = 1000 km, κ = 12 dB, B = 10 GBd, β = 0.05, EDFA on).
Run with GNU Octave 8.4, November 2026.

**Reproducibility.** Analytic numbers (noise budget, rates, sensitivity) are
deterministic and should match MATLAB to the digits shown. Monte Carlo numbers
depend on the random-number generator: MATLAB and Octave give different
sequences for the same seed, so expect differences of about ±0.1 dB in waveform
SNR. PAPR is the 10⁻⁴ point of the power CCDF over 2¹⁷ symbols, which varies
by only about ±0.1 dB between generators (the single largest peak used before
v2.5 varied by more than 1 dB).

## 1. Geometry and noise budget at y = A₀ (`RUN.budget`)

| Quantity | Value |
|---|---|
| w_z / σ_s / α | 4.870 m / 5.437 m / 0.2006 |
| A₀ | 4.132 × 10⁻⁴ |
| Received power (carrier / signal) | −9.52 dBm (−9.79 / −21.79 dBm) |
| ζ = 1 + 2/κ | 1.126 |
| SNR, shot only | 47.4 dB |
| SNR, thermal only | 48.9 dB |
| SNR, RIN only | 40.0 dB |
| SNR, signal/carrier–ASE only | 31.9 dB |
| SNR, ASE×ASE only | 74.8 dB |
| PAPR (10⁻⁴ CCDF), BPSK / OOK / QPSK | 7.5 / 7.4 / 6.4 dB |
| Heuristic reconstruction term, QPSK / BPSK (**not included**, `P.recon = false`) | 62.5 dB / 34.2 dB |
| **SNR total, pre-amp on (all modulations)** | **31.05 dB** |
| SNR total, photodiode only | 18.65 dB |

## 2. Rates and sensitivity

| Modulation | Ergodic, hard | Ergodic, soft | Sensitivity (BER 3.8×10⁻³) |
|---|---|---|---|
| BPSK | 6.84 Gb/s | 6.96 Gb/s | −30.4 dBm |
| OOK | 6.59 Gb/s | 6.72 Gb/s | −28.6 dBm |
| QPSK | **13.18 Gb/s** | **13.44 Gb/s** | −28.6 dBm |

- QPSK without pre-amplifier: 6.70 Gb/s (hard). BPSK: 3.80 Gb/s.
- 10 % outage rate: 7.6 × 10⁻⁵ Gb/s, effectively zero. The link is below the
  QPSK sensitivity 41 % of the time: (P_sens / P_rx,max)^α = (10^(−19.08/10))^0.2006.

| Detector radius w_d (`RUN.wd`) | 4 | 5 | 6 | 7 | 8 | 9 | 10 cm |
|---|---|---|---|---|---|---|---|
| QPSK, Gb/s | 11.46 | 12.19 | 12.74 | **13.18** | 13.53 | 13.83 | 14.09 |
| BPSK, Gb/s | 6.05 | 6.38 | 6.64 | 6.84 | 7.01 | 7.14 | 7.26 |

| C-band edge (`RUN.itu`) | w_z | A₀ | α | QPSK at 0.3 W | QPSK at 0.1 W |
|---|---|---|---|---|---|
| 1530 nm | 4.870 m | 4.13 × 10⁻⁴ | 0.201 | 13.18 Gb/s | 11.49 Gb/s |
| 1565 nm | 4.982 m | 3.95 × 10⁻⁴ | 0.210 | 13.45 Gb/s (+2.0 %) | 11.75 Gb/s (+2.3 %) |

## 3. Validation

**Channel Monte Carlo** (`RUN.chanMC`, 2 × 10⁵ Rayleigh samples):

| | MC | Analytic |
|---|---|---|
| E[y] | 6.927 × 10⁻⁵ | 6.904 × 10⁻⁵ |
| QPSK ergodic rate | 13.193 Gb/s | 13.176 Gb/s |
| BPSK ergodic rate | 6.846 Gb/s | 6.841 Gb/s |

**Waveform Monte Carlo** (`RUN.waveMC`, QPSK, 2¹³ symbols, 512 pilots):

| y/A₀ | 1 | 0.3 | 0.1 | 0.03 | 0.01 |
|---|---|---|---|---|---|
| Waveform SNR | 31.06 | 26.07 | 20.90 | 14.25 | 6.65 dB |
| Analytic g(y) | 31.05 | 26.05 | 20.91 | 14.33 | 7.05 dB |
| Waveform BER | 0 | 0 | 0 | 0 | 1.20 × 10⁻² |
| Analytic BER | — | — | — | 9.7 × 10⁻⁸ | 1.22 × 10⁻² |

**Noise terms tested one at a time** (waveform with only that noise on):

| Test | κ | Waveform | Analytic |
|---|---|---|---|
| Thermal only, no EDFA | 6 / 9 / 12 / 15 dB | 17.53 / 20.12 / 18.73 / 16.31 dB | 22.20 / 20.91 / 18.87 / 16.38 dB |
| ASE only, EDFA on | 9 / 12 / 15 dB | 34.32 / 31.91 / 29.19 dB | 34.35 / 31.85 / 29.11 dB |
| Noiseless, QPSK (reconstruction limit) | 6 / 9 / 12 dB | 44.8 / 78.2 / 100.9 dB | 10.5 / 27.7 / 62.0 dB (heuristic formula, single-peak PAPR 6.8 dB; +0.4 dB with the CCDF PAPR) |

So the ASE term is exact, ζ holds for κ ≥ 12 dB, and the heuristic reconstruction
formula is far too pessimistic. Since v2.5 it is off by default (`P.recon`).
Rates are unaffected (the term was negligible for QPSK at 12 dB and only lowered
the on-axis BPSK/OOK SNR from 31.05 to about 29 dB).

**Waveform SNR vs CSPR** (`RUN.cspr_wave`, QPSK):

| y/A₀ | Best κ | SNR at best κ | SNR at κ = 12 dB |
|---|---|---|---|
| 1 | 7 dB | 34.2 dB | 31.1 dB |
| 0.1 | 7 dB | 23.3 dB | 20.9 dB |
| 0.03 | 9 dB | 16.2 dB | 14.3 dB |

## 4. Constellation study at y = A₀ (`RUN.constellation`, QPSK)

7680 payload symbols (15 360 bits) after a 512-symbol pilot estimate.

| Configuration | SNR | EVM | BER | σ_φ | Power on diode |
|---|---|---|---|---|---|
| A2 Tx carrier | 31.06 dB | 2.8 % | 0 / 15 360 | 0.020 rad | 5.5 dBm (after EDFA) |
| B1 Rx LO (1.77 mW, ideal, no EDFA) | 37.30 dB | 1.4 % | 0 / 15 360 | 0.009 rad | 2.7 dBm |
| A1/B3 no carrier | −27.9 dB | > 100 % | 0.505 | 1.82 rad | 5.5 dBm |
| B2 carrier added after detection | −24.4 dB | > 100 % | 0.507 | 1.83 rad | 5.5 dBm |

B1 is 6.2 dB better than A2 with an ideal LO: the CSPR penalty of the transmitter
carrier, made visible. A realistic LO (offset, linewidth) reduces that margin;
set `MC.lo_df` and `MC.lo_lw` to study it.

## 5. CKK / GKK / MP-KK (`RUN.mpkk`, sensitivity study)

Imposed front-end distortion S_NL = S_ME = 10 dB (a₂ = 0.336, b = 0.207, tap
spacing 2 samples), QPSK, SNR in dB:

| y/A₀ | 1 | 0.5 | 0.2 | 0.1 | 0.05 | 0.03 | 0.02 | 0.01 |
|---|---|---|---|---|---|---|---|---|
| No distortion | 31.1 | 28.3 | 24.2 | 20.9 | 17.2 | 14.3 | 11.7 | 6.6 |
| CKK | 9.8 | 12.4 | 13.2 | 12.8 | 12.0 | 10.7 | 9.0 | 4.4 |
| CKK + post-KK EQ | 11.4 | 16.3 | 19.2 | 18.3 | 15.9 | 13.4 | 10.8 | 5.4 |
| GKK | 11.0 | 12.8 | 13.3 | 13.1 | 12.4 | 11.3 | 9.9 | 6.6 |
| MP-KK (K = 3, M = 7) | 16.4 | 24.1 | 24.0 | 20.9 | 17.3 | 14.4 | 11.9 | 7.3 |

BER is counted directly only where errors occur (y = 0.01·A₀ here); elsewhere
the validation rests on the measured SNR.

MP-KK vs K (M = 7): 11.1, 14.9, 16.4, 16.9, 17.2 dB for K = 1…5.
MP-KK vs M (K = 3): 11.0, 14.0, 15.7, 16.0, 16.4, 16.5, 16.5 dB for M = 1, 2, 3, 5, 7, 9, 13.

**Without the imposed distortion** (the actual link model) all receivers give the
same result at κ = 12 dB: 31.05 (CKK), 31.02 (CKK + EQ), 31.07 (GKK) and
31.05 dB (MP-KK) on-axis.

## 6. Residual Doppler and recovery (`RUN.doppler`, Section 11)

With a transmitter carrier the common Doppler shift cancels, but at v = 15 km/s
(±9.9 GHz carrier Doppler) the data sideband keeps a residual offset
f_sh·v/c ≈ 263 kHz and the symbol clock is scaled by v/c = 50 ppm. The receiver
estimates the offset from pilots (512 at the frame start, 128 at the end), takes
the clock correction from the same estimate (ε = Δf/f_sh), and tracks the
remaining phase blindly (4th power, 129-symbol window, leave-one-out). QPSK SNR
in dB:

| y/A₀ | 1 | 0.3 | 0.1 | 0.03 | 0.0123 (FEC threshold) | 0.01 |
|---|---|---|---|---|---|---|
| No Doppler | 31.06 | 26.07 | 20.90 | 14.25 | 8.28 | 6.65 |
| Doppler, no recovery | 1.99 | 1.98 | 1.95 | 1.77 | 1.12 | 0.78 |
| Doppler + recovery | 30.83 | 25.90 | 20.75 | 14.10 | 7.93 | 6.10 |
| Doppler (−v) + recovery | 30.97 | 25.95 | 20.76 | 14.08 | 8.10 | 6.33 |
| Estimated offset [kHz] | 266 | 268 | 272 | 284 | 317 | 338 |

Without recovery the residual Doppler destroys the link. With recovery the
penalty is ≤ 0.25 dB down to 0.03·A₀ and ≤ 0.4 dB at the FEC threshold (checked
over 6 random seeds). Not modelled: Tx/Rx oscillator mismatch, which needs
timing recovery in any receiver; Doppler rate (acceleration) within a frame.

## 7. Receiver-LO comparison with a realistic LO (`RUN.lo_compare`, Section 12)

B1: data at full P_T, LO 1.77 mW before the photodiode, no EDFA. The realistic LO
has a 100 MHz residual offset after locking and 200 kHz or 2 MHz combined
linewidth. B1 gets the same frequency/phase recovery, with the tracking window
chosen per case to give B1 its best result (optimistic for B1). QPSK SNR in dB:

| y/A₀ | 1 | 0.3 | 0.1 | 0.03 |
|---|---|---|---|---|
| A2 Tx carrier (Doppler + recovery) | 30.83 | 25.90 | 20.75 | 14.10 |
| B1 ideal LO | 37.30 | 32.74 | 28.15 | 22.98 |
| B1, 100 MHz offset, 200 kHz linewidth | 35.04 | 31.70 | 27.63 | 22.76 |
| B1, 100 MHz offset, 2 MHz linewidth | 31.36 | 29.55 | 26.40 | 22.14 |

A realistic LO reduces B1's on-axis advantage from 6.5 dB to 4.2 dB (200 kHz) or
0.5 dB (2 MHz), but at 0.1·A₀ and below B1 stays 5.7–8.7 dB better: its fixed LO power keeps the
beat signal above thermal noise (SNR ∝ y), while A2's carrier fades with the
signal (SNR ∝ y² when thermal-limited). The case for A2 therefore rests on the
hardware it avoids (LO laser, ±9.9 GHz tracking, frequency locking, polarisation
control), not on SNR.

## 8. Realistic pointing: power and range for a target outage (`RUN.pointing`, Section 13)

QPSK, z = 1000 km, κ = 12 dB. P(P_rx < P_sens) = (P_sens/(h_PL P_T A₀))^α, so
P_T,min = P_sens/(h_PL A₀ ε^{1/α}). Sensitivity −28.6 dBm with EDFA, −14.7 dBm
without.

| σ_s0 | α | Ergodic rate at 0.3 W | P_T for 1 % outage (EDFA / none) | P_T for 0.1 % outage (EDFA / none) |
|---|---|---|---|---|
| 0.25 m | 12.8 | 20.00 Gb/s | 7.2 / 21.2 dBm | 8.0 / 22.0 dBm |
| 0.5 m | 3.21 | 20.00 Gb/s | 11.9 / 25.9 dBm | 15.0 / 29.0 dBm |
| 0.75 m | 1.43 | 19.99 Gb/s | 19.7 / 33.6 dBm | 26.7 / 40.7 dBm |
| 1 m | 0.80 | 19.70 Gb/s | 30.6 / 44.6 dBm | 43.1 / 57.0 dBm |
| 2 m (stress case) | 0.20 | 13.18 Gb/s | 105 / 119 dBm (infeasible) | 155 / 169 dBm (infeasible) |

P_T = 0.3 W is 24.8 dBm. With realistic jitter (σ_s0 ≤ 0.75 m) QPSK saturates at
20 Gb/s, and the useful figure of merit is the power margin: about 13 dB for 1 %
outage at σ_s0 = 0.5 m. The pointing fade margin for 1 % outage is
ε^(−1/α): 6.2 dB at σ_s0 = 0.5 m and 1.6 dB at 0.25 m. Below about 0.5 m the rest of
the requirement is set by receiver sensitivity, which is where the CSPR penalty of
the transmitter carrier matters.

Feasible range window at P_T = 0.3 W (EDFA):

| σ_s0 | 1 % outage | 0.1 % outage |
|---|---|---|
| 0.25 m | 100–3000 km | 100–2750 km |
| 0.5 m | 200–2050 km | 300–1750 km |
| 0.75 m | 400–1400 km | none |
| ≥ 1 m | none | none |

The window has a lower edge because in the jitter model of Shang et al. σ_s is a
length that grows only slowly with range, while the beam radius grows in
proportion to range; at short range the beam is narrower than the jitter. With
angular jitter (σ_s = z·σ_θ) α would not depend on range (see
`open_questions.md`).

## 9. BER near the FEC threshold (`RUN.berfec`, Section 14)

QPSK, 8 runs × 2¹⁷ symbols (2.09 × 10⁶ bits per point), pointing states chosen
where the analytic BER equals the target:

| Target (analytic) | y/A₀ | Counted BER (95 % interval) | Errors |
|---|---|---|---|
| 3.8 × 10⁻³ (HD-FEC) | 0.0123 | 3.76 × 10⁻³ ± 0.08 × 10⁻³ | 7849 |
| 1.0 × 10⁻⁴ | 0.0188 | 0.98 × 10⁻⁴ ± 0.13 × 10⁻⁴ | 205 |

The BER, not only the SNR, matches the analytic model where it matters for FEC.
