# The KK-OISL model

Summary of the model implemented in `matlab/oisl_kk_txcarrier.m`. Equation numbers
refer to Chapter 3 of the report. Numbers are for the default parameters.

## 1. Transmitter and channel

Carrier and SSB data come from one laser and one IQ modulator:

$$E_{Tx}(t) = \Big[\sqrt{P_c} + \sqrt{P_s}\,s(t)\,e^{i2\pi f_{sh}t}\Big]e^{i\theta_L(t)}, \qquad f_{sh} = (1+\beta)B/2$$

$$P_c = \frac{\kappa}{1+\kappa}P_T, \qquad P_s = \frac{1}{1+\kappa}P_T, \qquad \kappa = \text{CSPR}$$

Both parts see the same channel gain $h = h_{PL}\,y$, where $y \in (0, A_0]$ is
the pointing fraction, so the CSPR at the receiver equals κ for every pointing
state. KK needs the field to be minimum-phase; a sufficient condition is
$\kappa > \mathrm{PAPR}_s$.

## 2. Receiver chain (waveform simulation, `kk_link`)

1. Optional EDFA pre-amplifier: gain $G$ = 15 dB, noise figure 5 dB, ASE in both
   polarisations, optical band-pass filter $B_f$ = 30 GHz.
2. Photodiode: $I = \eta(|E_x|^2 + |E_y|^2)$ plus shot, thermal and RIN noise;
   electrical bandwidth $B_o = (1+\beta)B$.
3. KK DSP: $\phi = \mathcal{H}\{\tfrac12\ln I_n\}$,
   $\hat E = \sqrt{I_n}\,e^{i\phi}$, remove the carrier (mean), shift down by
   $f_{sh}$.
4. RRC matched filter, one sample per symbol.
5. Pilot-aided one-tap gain/phase estimate from the first 512 symbols. SNR, BER,
   EVM and phase error are measured on the remaining 7680 payload symbols.

## 3. Analytic SNR (`rx_snr`)

With $P_{c,d}$, $P_{s,d}$ the carrier and signal powers at the photodiode,
$S = 2\eta^2 P_{c,d}P_{s,d}$ and the KK noise factor $\zeta = 1 + 2/\kappa$:

| Term | SNR contribution |
|---|---|
| Shot, thermal, RIN, ASE×ASE | $g_j = S / (\zeta N_j B)$ |
| Signal/carrier–ASE | $g_{ASE} = P_{s,d} / [(1+\zeta) S_{ASE} B]$ |
| Reconstruction penalty (heuristic) | $1/g_{rec} = \mathrm{PAPR}_s \cdot P_{viol}(\kappa)$ |

$$g(y) = \Big(\sum_j 1/g_j\Big)^{-1}$$

Validated term by term against the waveform simulation (see `results.md`):
the ASE term and ζ hold within 0.15 dB for κ ≥ 12 dB. The reconstruction penalty
is far more pessimistic than the waveform at low κ and should not be used to pick
the CSPR; the waveform sweep (`RUN.cspr_wave`) is the reference for that.

## 4. Pointing fading and rates

$$w_z = z\tan\frac{\lambda}{\pi w_0},\quad A_0 = \frac{2w_d^2}{w_z^2},\quad \sigma_s = \sigma_{s0}e^{k_0 z/d_0},\quad \alpha = \frac{w_z^2}{4\sigma_s^2}$$

$$f_Y(y) = \frac{\alpha}{A_0^\alpha}y^{\alpha-1}, \qquad F_Y(y) = (y/A_0)^\alpha$$

- Ergodic rate: $R = B\int_0^{A_0} m(g(y)) f_Y(y)\,dy$, with $m$ the hard- or
  soft-decision mutual information of BPSK, OOK or QPSK.
- ε-outage rate: $R_\epsilon = B\,m\big(g(A_0\epsilon^{1/\alpha})\big)$.
- Sensitivity: received power at which BER = 3.8×10⁻³ (7 % HD-FEC).

## 5. Carrier/LO placements (constellation study, Section 8)

| Label | Config | What reaches the photodiode |
|---|---|---|
| A2 | `tx_carrier` | carrier + data, both faded by y (the project's design) |
| B1 | `rx_lo` | data at full P_T + local LO before the diode; no EDFA; ideal LO by default (`MC.lo_df`, `MC.lo_lw` add offset and linewidth) |
| A1/B3 | `no_carrier` | data only |
| B2 | `post_det_lo` | data only; the carrier level is added digitally after detection |

## 6. Optional: CKK / GKK / MP-KK (Section 10, `RUN.mpkk`)

Following Wen et al. (Opt. Express 2023), an **imposed** front-end distortion
acts on the photocurrent: nonlinearity $v = u + a_2u^2$ followed by memory
$w[n] = v[n] + b\,v[n-d] + b\,v[n-2d]$, with $a_2$, $b$ set from $S_{NL}$, $S_{ME}$
in dB. Receivers compared: CKK, CKK + post-KK equaliser, GKK (memoryless
polynomial before KK) and MP-KK (memory polynomial before KK, trained on the
pilots). This is a sensitivity study: the link model has no nonlinearity or memory
of its own.
