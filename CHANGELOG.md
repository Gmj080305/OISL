# Changelog

## November 2026: repository restructure

- New layout: `matlab/`, `figures/`, `report/`, `docs/`, `references/`, `slides/`.
- The README now describes the KK-OISL project; the original notes moved to
  `docs/background.md`, with short corrections added at the top.
- The base-paper PDF was replaced by links in `references/README.md` (publisher
  copyright). The NotebookLM slide decks moved to `slides/`.
- Added `docs/kk_oisl_model.md`, `docs/results.md`, `docs/open_questions.md`
  and `report/REVISIONS.md`.

### matlab/oisl_kk_txcarrier.m

**v2.4 (Nov 2026)**
- Every figure is saved when `PLOT.save = true` (PNG + PDF, in `figures/`).
- `RUN.wd` (Fig 3.8, rate vs detector radius) and `RUN.itu` (Fig 3.9, rate vs
  P_T at 1530 / 1565 nm) added, so all report figures come from MATLAB.
- `RUN.mpkk` is off by default.
- Bug fix: in Fig 3.4 the reconstruction term was a scalar and was drawn as 200
  separate lines, filling the legend. It now has the size of the power vector
  (no change to any computed number).

**v2.3 (Oct 2026)**
- The constellation study is drawn in one black-background window: per
  configuration the constellation (scatter, coloured by sent symbol), a density map
  and a phase-error histogram, plus the A2 overlay across pointing states.
- Optional `scatterplot()` overlay (Communications Toolbox) in its own window.

**v2.2 (Oct 2026)**
- Section 10: CKK, CKK + post-KK equaliser, GKK and MP-KK receivers with an
  optional imposed front-end distortion (Wen et al., Opt. Express 2023).

**v2.1 (Oct 2026)**
- Bug fix: `P.c` is now defined before `P.f = P.c/1530e-9` (the script stopped at
  that line).
- The waveform Monte Carlo and the constellation study share one waveform core
  (`kk_link`), so their numbers agree (A2: 31.1 dB, the analytic value).
- Pilot-aided gain/phase estimate (512 symbols); metrics on payload only.
- B1 (receiver LO) now matches its description: data at full P_T, LO before the
  diode, no EDFA, optional LO offset and linewidth.
- B2 is now "carrier added after detection", which is physically different from
  A1/B3 (previously the two gave identical numbers).
- `RUN.cspr_wave`: waveform SNR vs CSPR.
- MATLAB / Octave compatibility helpers; no toolbox needed.

### matlab/oisl_kk_channel.m

- Unknown parameter names, an odd number of arguments and unsupported `mod`
  values raise an error (a typo such as `'Pt'` used to be ignored silently).
- `config` other than `'tx_carrier'` raises an error: only that case has an
  analytic model. The header no longer claims four configurations.
- The PAPR calculation no longer resets the caller's random-number state.
- Header lists every output field.

### New files

- `matlab/run_in_octave.m`: runs the main script in GNU Octave.
- `matlab/tests/check_oisl_kk.m`: 16 regression checks for `oisl_kk_channel.m`.

## June 2026

- Initial repository: README with OISL system-model notes and a summary of
  Shang et al. (2025); base paper PDF; NotebookLM slides.
