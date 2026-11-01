# Reference results

Numbers produced by the code in `matlab/` (v2.4) with the default parameters
(1530 nm, P_T = 0.3 W, z = 1000 km, κ = 12 dB, B = 10 GBd, β = 0.05, EDFA on).
Run with GNU Octave 8.4, November 2026.

**Reproducibility.** Analytic numbers (noise budget, rates, sensitivity) are
deterministic and should match MATLAB to the digits shown. Monte Carlo numbers
depend on the random-number generator: MATLAB and Octave give different
sequences for the same seed, so expect differences of about ±0.1 dB in waveform
SNR. The measured PAPR is the single largest peak of one random waveform and can
differ by more than 1 dB between generators (Octave: QPSK 6.8 dB; the earlier
Python runs: 8.5 dB). At κ = 12 dB this has no visible effect on QPSK results.

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
| SNR, reconstruction term, QPSK / BPSK | 62.0 dB / 33.7 dB (PAPR-dependent) |
| **SNR total, QPSK, pre-amp on** | **31.05 dB** |
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
| Noiseless, QPSK (reconstruction limit) | 6 / 9 / 12 dB | 44.8 / 78.2 / 100.9 dB | 10.5 / 27.7 / 62.0 dB (formula) |

So the ASE term is exact, ζ holds for κ ≥ 12 dB, and the reconstruction formula
is far too pessimistic.

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

MP-KK vs K (M = 7): 11.1, 14.9, 16.4, 16.9, 17.2 dB for K = 1…5.
MP-KK vs M (K = 3): 11.0, 14.0, 15.7, 16.0, 16.4, 16.5, 16.5 dB for M = 1, 2, 3, 5, 7, 9, 13.

**Without the imposed distortion** (the actual link model) all receivers give the
same result at κ = 12 dB: 31.05 (CKK), 31.02 (CKK + EQ), 31.07 (GKK) and
31.05 dB (MP-KK) on-axis.
