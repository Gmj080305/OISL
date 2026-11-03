%% ========================================================================
%  KK-OISL CHANNEL SIMULATION & RATE ANALYSIS  (version 2)
%  Self-coherent KK receiver: the CW carrier is sent from the TRANSMITTER
%  together with the single-sideband (SSB) data signal.
%  Modulations: BPSK / OOK / QPSK
%  ========================================================================
%  Model (Chapter 3, eqs. 3.1-3.22):
%   Tx field     E_T = sqrt(Pc) + sqrt(Ps)*s(t)*exp(i*2*pi*fsh*t),  fsh = Bo/2
%                Pc = k/(1+k)*PT,  Ps = PT/(1+k),   k = CSPR (linear)
%   Rx power     Prx = hPL*y*PT  (carrier and signal scaled by the same y)
%   Pre-amp      optional EDFA, fixed gain G, noise figure F (G = 1 -> off)
%                ASE PSD per polarisation  S = h*nu*(F*G-1)/2
%   Photocurrent noise (one-sided PSDs, A^2/Hz):
%                shot 2*q*Ibar, thermal i_th^2, RIN*(eta*(Pc+Ps))^2,
%                ASE x ASE 4*eta^2*S^2*(Bf-B/2)
%   SNR mapping  g_pc  = 2*eta^2*Pc_d*Ps_d / (zeta*N_i*B),  zeta = 1+2/k
%                g_ase = Ps_d / ((1+zeta)*S*B)
%                g = 1/(sum of 1/g_j)
%   Rates        hard: m = bits*(1-H2(BER)); soft: Gaussian-input-free MI (GMI)
%                ergodic R = B*int m(g(y)) f_Y(y) dy, f_Y = a*A0^-a*y^(a-1)
%                outage  R_eps = B*m(g(A0*eps^(1/a)))
%   Monte Carlo  (i) channel: Rayleigh pointing error -> y -> rate average
%   Monte Carlo  (i) channel: Rayleigh pointing error -> y -> rate average
%                (ii) waveform: full KK link at chosen y -> SNR and BER (pilot-aided)
%                (iii) constellations: 4 carrier/LO placements, same waveform core
%   v2.1 (Oct 2026): P.c defined before use; constellation study rewritten
%                (shared kk_link core, pilot-aided metrics, B1 = LO before the
%                diode with no EDFA, B2 = carrier added after detection, colour /
%                density / pointing-overlay plots); optional waveform CSPR sweep.
%   v2.6 (Nov 2026): Sections 11-14: residual Doppler and recovery DSP, fair
%                receiver-LO comparison, realistic-pointing power/range study,
%                BER near the FEC threshold.
%   v2.5 (Nov 2026): heuristic reconstruction penalty off by default (P.recon);
%                PAPR = 10^-4 CCDF point; waveform core moved to function files
%                (kk_link.m, rx_metrics.m, const_map.m, rrc.m, fe_calibrate.m).
%   v2.4 (Nov 2026): every figure saved with PLOT.save; Fig 3.4 recon term fixed
%                (was a scalar, drawn as 200 separate lines); RUN.wd (Fig 3.8) and
%                RUN.itu (Fig 3.9) added; RUN.mpkk off by default.
%   v2.3 (Oct 2026): Section 8 figures merged into one black-background window;
%                point clouds drawn with scatter() + hold on (scatterplot() optional).
%   v2.2 (Oct 2026): Section 10, CKK / CKK+EQ / GKK / MP-KK receivers with an
%                optional front-end distortion model (Wen et al., Opt. Express 2023).
clear; clc; close all;

%% ===================== 1. USER PARAMETERS =================================
% ---- constants (defined first: P.f below uses P.c) ----
P.c = 299792458; P.q = 1.602176634e-19; P.h = 6.62607015e-34;
% ---- geometry / channel (same as Chapter 2) ----
P.w0       = 0.1;        % [m]    beam waist at Tx
P.wd       = 0.07;       % [m]    detector (aperture) radius (20-50% < beam waist)
P.f        = P.c/1530e-9; % [Hz]  default: ITU-T C-band short edge (1530 nm)
P.hPL      = 0.9;        % [-]    terminal / path efficiency
P.PT       = 0.3;        % [W]    default total transmit power (carrier + signal)
P.z        = 1000e3;     % [m]    default link distance
P.sigma_s0 = 2;          % [m]    jitter model sigma_s(z) = sigma_s0*exp(k0*z/d0)
P.k0       = 0.1;
P.d0       = 100e3;      % [m]
% ---- signal ----
P.B        = 10e9;       % [Hz]   symbol rate (= signal bandwidth for Nyquist pulses)
P.CSPR_dB  = 12;         % [dB]   carrier-to-signal power ratio at Tx
P.beta     = 0.05;       % [-]    RRC roll-off factor (used for PAPR computation)
% ---- receiver ----
P.eta      = 0.5;        % [A/W]  photodiode responsivity
P.i_th     = 20e-12;     % [A/sqrt(Hz)] thermal noise current density
P.RIN_dB   = -150;       % [dB/Hz] transmit-laser RIN
P.pre_on   = true;       % true = EDFA pre-amplifier before the photodiode
P.G_dB     = 15;         % [dB]   pre-amp gain (fixed)
P.NF_dB    = 5;          % [dB]   pre-amp noise figure
P.Bf       = 30e9;       % [Hz]   optical band-pass filter after the pre-amp
P.pth      = 0;          % [W]    optional tracking threshold (0 = none)
P.recon    = false;      % true: add the heuristic KK reconstruction penalty
                         %       PAPR*P_viol(kappa). It is 35-50 dB more pessimistic than
                         %       the waveform simulation (docs/results.md), so it is off by
                         %       default and only kept for comparison.

MODS = {'BPSK','OOK','QPSK'};     % modulations to evaluate

% ---- sweep ranges (edit freely) ----
SW.f_vec     = linspace(191e12, 196e12, 25);   % [Hz]  Fig 3.1 / 3.2 (ITU-T C-band)
SW.PT_vec    = 0.10:0.02:0.30;                 % [W]   Fig 3.1 / 3.2 curves
SW.SNR_dB    = -5:0.25:30;                     % [dB]  Fig 3.3
SW.yrel      = logspace(-4, 0, 200);           % [-]   Fig 3.4 (y/A0)
SW.CSPR_dB   = 0:0.5:20;                       % [dB]  CSPR sweep (wide range to show bell curve)
SW.sig0_vec  = linspace(0.2, 4, 30);           % [m]   optional jitter sweep
SW.z_vec     = [500e3 1000e3 1500e3];          % [m]   optional jitter sweep curves
SW.eps_vec   = logspace(-3, -0.3, 40);         % [-]   optional outage sweep
SW.wd_vec    = 0.04:0.01:0.10;                 % [m]   detector radius sweep (Fig 3.8)
SW.v_max     = 15e3;                           % [m/s] relative radial velocity (+-9.9 GHz Doppler)
SW.trk_win_vec = [3 5 9 17 33 65];             % [sym] phase-tracking windows tried for B1
SW.sig0_real = [0.25 0.5 0.75 1 1.5 2];        % [m]   jitter values for the pointing study
SW.avail_eps = [1e-2 1e-3];                    % [-]   outage targets (99 %, 99.9 % availability)
SW.z_range   = (100:50:5000)*1e3;              % [m]   range grid for the maximum-range search
SW.ber_targets = [3.8e-3 1e-4];                % [-]   analytic BER points for the BER count
SW.ber_runs  = 8;  SW.ber_log2N = 17;          % 8 x 2^17 symbols (about 1e6) per point

% ---- Monte Carlo settings ----
MC.Nchan     = 2e5;               % channel MC samples
MC.wave_mod  = 'QPSK';            % modulation for waveform MC
MC.wave_yrel = [1 0.3 0.1 0.03 0.01];  % y/A0 points for waveform MC
MC.Nsym      = 2^13;              % symbols per waveform run
MC.Nos       = 8;                 % samples per symbol
MC.beta      = P.beta;            % RRC roll-off (from signal params)
MC.seed      = 1;
MC.Npilot    = 512;               % known pilot symbols (gain/phase estimate, compensator training)
MC.guard     = 16;                % pilot symbols at each block edge left out of compensator training
% ---- B1 (Rx optical LO) settings for the constellation study ----
MC.lo_kappa  = 10^(12/10);        % LO-to-signal power ratio at y = MC.lo_yref (set below to A0)
MC.lo_df     = 0;                 % residual LO frequency offset after locking [Hz]  (0 = ideal LO)
MC.lo_lw     = 0;                 % combined Tx + LO laser linewidth [Hz]           (0 = ideal LO)
MC.lo_real_df = 100e6;            % realistic B1 LO: residual offset after locking [Hz] (Section 12)
MC.lo_real_lw = [200e3 2e6];      % realistic B1 LO: combined linewidths to compare [Hz]
MC.Npost     = 128;               % known symbols at the frame end (used only with recovery)
MC.trk_win   = 129;               % phase-tracking window [symbols] (A2; B1 window is optimised)
SW.CSPR_wave_dB = 3:1:15;         % [dB]  waveform CSPR sweep (RUN.cspr_wave)
SW.cspr_yrel    = [1 0.1 0.03];   % [-]   pointing states for the waveform CSPR sweep

% ---- figure output ----
PLOT.use_scatterplot = true;      % also draw the pointing overlay with Communications Toolbox
                                  % scatterplot() + hold on (own window; skipped if not installed)
PLOT.save = false;                % true: save every figure as PNG + PDF in PLOT.dir
PLOT.dir  = fullfile('..', 'figures');   % repo-level figures/ folder

% ---- compute signal PAPR for each modulation (needed for KK penalty) ----
for mm = 1:numel(MODS)
    P.PAPR.(MODS{mm}) = signal_papr(MODS{mm}, P.B, P.beta);
end
fprintf('Signal PAPR (10^-4 CCDF): BPSK %.1f dB, OOK %.1f dB, QPSK %.1f dB\n', ...
    10*log10(P.PAPR.BPSK), 10*log10(P.PAPR.OOK), 10*log10(P.PAPR.QPSK));

% ---- which sections to run ----
RUN.budget   = true;   % Table 3.2 noise budget printout
RUN.freq     = true;   % Fig 3.1 / 3.2  rate vs frequency over PT
RUN.snr      = true;   % Fig 3.3        rate vs SNR (hard + soft)
RUN.ybudget  = true;   % Fig 3.4        SNR contributions vs y
RUN.chanMC   = true;   % Fig 3.5a       channel Monte Carlo
RUN.waveMC   = true;   % Fig 3.5b       waveform KK Monte Carlo
RUN.cspr     = false;  % optional: rate vs CSPR
RUN.jitter   = false;  % optional: rate vs sigma_s0 over z
RUN.outage   = false;  % optional: outage rate vs eps
RUN.constellation = true;  % Fig 3.6 / 3.7: constellation study (4 carrier/LO configs)
RUN.cspr_wave = false;  % optional: waveform SNR vs CSPR (true optimum CSPR, slower)
RUN.wd       = false;  % optional: Fig 3.8  rate vs detector radius
RUN.itu      = false;  % optional: Fig 3.9  rate vs PT at 1530 / 1565 nm
RUN.mpkk     = false;  % optional: Fig 3.10-3.13 CKK / CKK+EQ / GKK / MP-KK with an
                       %           IMPOSED front-end distortion (sensitivity study only)
RUN.doppler  = false;  % optional: residual Doppler with and without recovery (Section 11)
RUN.lo_compare = false; % optional: A2 vs receiver LO with realistic offset/linewidth (Section 12)
RUN.pointing = false;  % optional: required P_T and range for a target outage (Section 13)
RUN.berfec   = false;  % optional: counted BER near the FEC threshold, ~1e6 symbols (Section 14)
RUN.csv      = false;  % optional: CSV of rate vs f vs PT
csv_file     = 'rate_vs_freq_PT_v2.csv';

%% ===================== 2. NOISE BUDGET (Table 3.2) ========================
GQ = geom(P, P.f, P.z);
MC.lo_yref = GQ.A0;               % B1 LO power referenced to the on-axis received power
if RUN.budget
    Prx = P.hPL*P.PT*GQ.A0;
    [g, T] = rx_snr(P, P.f, Prx);
    k = 10^(P.CSPR_dB/10);
    fprintf('\n--- Noise budget: PT=%.2f W, CSPR=%.1f dB, %.0f THz, %.0f km, y=A0=%.3e ---\n', ...
        P.PT, P.CSPR_dB, P.f/1e12, P.z/1e3, GQ.A0);
    fprintf('  w_z = %.3f m, sigma_s = %.3f m, alpha = %.4f\n', GQ.wz, GQ.sig, GQ.alpha);
    fprintf('  Received power         : %.2f dBm (carrier %.2f dBm, signal %.2f dBm)\n', ...
        10*log10(Prx/1e-3), 10*log10(k/(1+k)*Prx/1e-3), 10*log10(Prx/(1+k)/1e-3));
    fprintf('  zeta (KK noise factor) : %.3f\n', 1+2/k);
    fn = fieldnames(T);
    for i = 1:numel(fn)
        fprintf('  SNR (%-8s only)    : %.1f dB\n', fn{i}, -10*log10(T.(fn{i})));
    end
    fprintf('  SNR (total)            : %.1f dB\n', 10*log10(g));
    for m = 1:numel(MODS)
        gm = rx_snr(P, P.f, Prx, MODS{m});
        fprintf('  %-4s PAPR = %.1f dB, heuristic recon SNR = %.1f dB (%s), total SNR = %.1f dB\n', ...
            MODS{m}, 10*log10(P.PAPR.(MODS{m})), -10*log10(recon_term(P, MODS{m})), ...
            onoff_str(P.recon, 'included', 'not included'), 10*log10(gm));
        sp = sensitivity(P,P.f,MODS{m});
        if isinf(sp), sens_str = 'inf (KK fails)';
        else, sens_str = sprintf('%.1f dBm', 10*log10(sp/1e-3)); end
        fprintf('  %-4s ergodic rate      : %.2f Gb/s (hard), %.2f Gb/s (soft); sensitivity %s\n', ...
            MODS{m}, R_erg(P,P.f,P.PT,MODS{m},'hard')/1e9, R_erg(P,P.f,P.PT,MODS{m},'soft')/1e9, sens_str);
    end
end

%% ===================== 3. Fig 3.1 / 3.2 : RATE vs FREQUENCY =================
if RUN.freq
    for m = [1 3]   % BPSK and QPSK as in the report (change to 1:numel(MODS) for all)
        R = zeros(numel(SW.PT_vec), numel(SW.f_vec));
        for i = 1:numel(SW.PT_vec)
            for j = 1:numel(SW.f_vec)
                R(i,j) = R_erg(P, SW.f_vec(j), SW.PT_vec(i), MODS{m}, 'hard');
            end
        end
        figure('Name',['Rate vs f ' MODS{m}]); hold on; grid on;
        cl = parula_compat(); cl = cl(round(linspace(1, 220, numel(SW.PT_vec))), :);  % ordered by P_T
        for i = 1:numel(SW.PT_vec)
            plot(SW.f_vec/1e12, R(i,:)/1e9, 'LineWidth',1.5, 'Color',cl(i,:), ...
                'DisplayName', sprintf('P_T = %.2f W', SW.PT_vec(i)));
        end
        xlabel('Laser Frequency f [THz]'); ylabel(['Ergodic Rate R_{KK-' MODS{m} '} [Gbps]']);
        title(['KK-' MODS{m} ' Rate vs Laser Frequency over Transmit Power P_T']);
        legend('Location','eastoutside');
        save_fig(PLOT, sprintf('fig3%d_rate_vs_freq_%s', 1 + (m == 3), MODS{m}));
    end
end

%% ===================== 4. Fig 3.3 : RATE vs SNR (AWGN) =====================
if RUN.snr
    gv = 10.^(SW.SNR_dB/10);
    figure('Name','Rate vs SNR'); hold on; grid on; cl = lines(numel(MODS));
    for m = 1:numel(MODS)
        plot(SW.SNR_dB, P.B*mi(MODS{m}, gv, 'hard')/1e9, '-',  'LineWidth',1.6, 'Color',cl(m,:), 'DisplayName',[MODS{m} ' hard']);
        plot(SW.SNR_dB, P.B*mi(MODS{m}, gv, 'soft')/1e9, '--', 'LineWidth',1.6, 'Color',cl(m,:), 'DisplayName',[MODS{m} ' soft']);
    end
    xlabel('SNR g = E_s/N_0 [dB]'); ylabel('Rate B\cdotm(g) [Gbps]');
    title('KK Rate vs SNR (hard and soft decision)'); legend('Location','southeast');
    save_fig(PLOT, 'fig33_rate_vs_snr');
end

%% ===================== 5. Fig 3.4 : SNR CONTRIBUTIONS vs y ==================
if RUN.ybudget
    y = GQ.A0*SW.yrel;
    [g, T] = rx_snr(P, P.f, P.hPL*P.PT*y, 'QPSK');    % recon term shown only if P.recon
    figure('Name','SNR budget vs y'); hold on; grid on;
    fn = fieldnames(T);
    for i = 1:numel(fn)
        lbl = fn{i}; if strcmp(lbl,'recon'), lbl = 'recon (QPSK)'; end
        semilogx(SW.yrel, -10*log10(T.(fn{i})), 'LineWidth',1.2, 'DisplayName', lbl);
    end
    semilogx(SW.yrel, 10*log10(g), 'k', 'LineWidth',2, 'DisplayName','total g(y)');
    set(gca,'XScale','log'); ylim([-20 70]);
    xlabel('y / A_0'); ylabel('SNR contribution [dB]');
    title(sprintf('SNR budget (QPSK), CSPR = %.0f dB, pre-amp %s', P.CSPR_dB, onoff(P.pre_on)));
    legend('Location','southeast');
    save_fig(PLOT, 'fig34_snr_budget_vs_y');
end

%% ===================== 6. Fig 3.5a : CHANNEL MONTE CARLO ====================
if RUN.chanMC
    rng(MC.seed);
    r  = GQ.sig*sqrt(-2*log(rand(MC.Nchan,1)));      % Rayleigh radial error
    ys = GQ.A0*exp(-2*r.^2/GQ.wz^2);                  % collected fraction (eq. 3.17 model)
    yq = sort(ys); Fe = (1:MC.Nchan)'/MC.Nchan;
    figure('Name','Channel MC');
    loglog(yq/GQ.A0, Fe, 'b', 'LineWidth',1.5); hold on; grid on;
    yy = logspace(-6,0,100); loglog(yy, yy.^GQ.alpha, 'r--', 'LineWidth',1.5);
    xlabel('y / A_0'); ylabel('F_Y(y)'); legend('Monte Carlo','(y/A_0)^\alpha','Location','southeast');
    title('Pointing-loss CDF: Monte Carlo vs analytic');
    save_fig(PLOT, 'fig35a_channel_mc_cdf');
    fprintf('\nChannel MC: E[y] = %.4e (MC) vs %.4e (analytic alpha*A0/(1+alpha))\n', ...
        mean(ys), GQ.alpha*GQ.A0/(1+GQ.alpha));
    for m = 1:numel(MODS)
        Rmc = P.B*mean(mi(MODS{m}, rx_snr(P,P.f,P.hPL*P.PT*ys,MODS{m}), 'hard'));
        fprintf('  %-4s ergodic rate: MC %.3f Gb/s vs integral %.3f Gb/s\n', MODS{m}, Rmc/1e9, ...
            R_erg(P,P.f,P.PT,MODS{m},'hard')/1e9);
    end
end

%% ===================== 7. Fig 3.5b : WAVEFORM KK MONTE CARLO ================
if RUN.waveMC
    nY = numel(MC.wave_yrel); snr_mc = zeros(1,nY); ber_mc = zeros(1,nY);
    g_an = zeros(1,nY); ber_an = zeros(1,nY);
    for i = 1:nY
        y = GQ.A0*MC.wave_yrel(i);
        o = kk_wave_mc(P, P.f, P.PT, y, MC.wave_mod, MC);
        snr_mc(i) = o.snr; ber_mc(i) = o.ber;
        g_an(i) = rx_snr(P, P.f, P.hPL*P.PT*y, MC.wave_mod); ber_an(i) = ber_theory(MC.wave_mod, g_an(i));
        fprintf('Wave MC %s y/A0=%.3g: SNR MC %.2f dB / analytic %.2f dB, BER MC %.2e / analytic %.2e\n', ...
            MC.wave_mod, MC.wave_yrel(i), 10*log10(o.snr), 10*log10(g_an(i)), o.ber, ber_an(i));
    end
    figure('Name','Waveform MC');
    subplot(1,2,1); semilogx(MC.wave_yrel, 10*log10(g_an),'r-o', MC.wave_yrel, 10*log10(snr_mc),'bs','LineWidth',1.4);
    grid on; xlabel('y / A_0'); ylabel('SNR [dB]'); legend('analytic g(y)','waveform MC','Location','southeast');
    subplot(1,2,2); loglog(MC.wave_yrel, max(ber_an,1e-12),'r-o', MC.wave_yrel, max(ber_mc,1e-12),'bs','LineWidth',1.4);
    grid on; xlabel('y / A_0'); ylabel('BER'); legend('analytic','waveform MC (0 errors -> 1e-12)', 'Location','southwest');
    title_compat(['End-to-end KK waveform Monte Carlo, ' MC.wave_mod]);
    save_fig(PLOT, 'fig35b_waveform_mc');
end

%% ===================== 8. CONSTELLATION STUDY (Fig 3.6 / 3.7) ==============
%  All four configurations share one waveform core (kk_link) and one receiver
%  metric (rx_metrics): pilot-aided one-tap gain/phase from the first
%  MC.Npilot symbols, then SNR / BER / EVM / phase error on the payload only.
%  Plots use a unit-power view so every panel shares the same axes.
if RUN.constellation
    CF.configs = {'tx_carrier','rx_lo','no_carrier','post_det_lo'};
    CF.labels  = {'A2: Tx carrier (co-propagating)', ...
                  'B1: Rx optical LO (before diode)', ...
                  'A1/B3: no carrier (KK fails)', ...
                  'B2: carrier added after diode (digital)'};
    CF.mod     = 'QPSK';
    CF.y       = GQ.A0;                    % on-axis: isolates carrier placement from fading
    CF.lim     = 1.8;                      % axis limit of the unit-power view
    [pts_c, ~] = const_map(CF.mod);
    pts_n      = pts_c/sqrt(mean(abs(pts_c).^2));
    symcol     = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];

    fprintf('\n--- Constellation study at y = A0, %s, %d payload symbols, %d pilots ---\n', ...
        CF.mod, MC.Nsym-MC.Npilot, MC.Npilot);
    fprintf('%-12s | %8s | %8s | %-22s | %9s | %7s | %9s\n', ...
        'Config','SNR[dB]','EVM[%]','BER (errors/bits)','sigma_phi','Pd[dBm]','P_LO[mW]');
    fprintf('%s\n', repmat('-',1,94));
    res = cell(1,4);
    for ci = 1:4
        L = kk_link(P, P.f, P.PT, CF.y, CF.mod, MC, struct('config', CF.configs{ci}));
        res{ci} = rx_metrics(L, MC.Npilot);
        o = res{ci};
        fprintf('%-12s | %8.2f | %8.2f | %9.2e (%5d/%5d) | %9.4f | %7.1f | %9.2f\n', ...
            CF.configs{ci}, 10*log10(o.snr), 100*o.evm, o.ber, o.nerr, o.nbits, ...
            o.sigma_phi, L.Pd_dBm, L.P_LO*1e3);
    end
    fprintf('B1 note: LO set for CSPR = %.0f dB at y = A0, no EDFA, ideal LO (lo_df = %g Hz, lo_lw = %g Hz).\n', ...
        10*log10(MC.lo_kappa), MC.lo_df, MC.lo_lw);

    % ---- Fig 3.6: everything in ONE window, black background ----
    %  Columns 1-4: one carrier/LO configuration each
    %     row 1  constellation, coloured by the transmitted symbol
    %     row 2  density (2-D histogram, log scale)
    %     row 3  phase-error distribution (payload only)
    %  Columns 5-6 (rows 1-2): A2 at several pointing states overlaid with HOLD ON
    symcol = [0.30 0.75 1.00; 1.00 0.55 0.20; 0.55 0.95 0.35; 0.95 0.45 0.95];  % bright on black
    hf = figure('Name','Constellation study','Color','k','Position',[30 40 1750 980]);
    if isprop(hf, 'InvertHardcopy'), set(hf, 'InvertHardcopy', 'off'); end   % keep black when saving
    edges = linspace(-CF.lim, CF.lim, 121);
    cmap  = [0 0 0; hot(255)];                 % zero counts stay black
    for ci = 1:4
        o = res{ci};
        % row 1: constellation
        ax = subplot(3,6,ci); hold on;
        for j = 1:numel(pts_c)
            sel = (o.idx == j);
            scatter(real(o.sym_plot(sel)), imag(o.sym_plot(sel)), 3, symcol(j,:), 'filled');
        end
        scatter(real(pts_n), imag(pts_n), 60, [1 1 1], '+', 'LineWidth', 1);
        axis equal; xlim(CF.lim*[-1 1]); ylim(CF.lim*[-1 1]);
        xlabel('In-phase'); ylabel('Quadrature');
        title([{CF.labels{ci}}, metrics_str(o, CF.mod)], 'FontSize', 8);
        dark_axes(ax);
        % row 2: density
        ax = subplot(3,6,6+ci);
        H2 = hist2d_counts(real(o.sym_plot), imag(o.sym_plot), edges);
        imagesc(edges, edges, log10(1 + H2)); axis xy; axis equal tight; hold on;
        plot(real(pts_n), imag(pts_n), 'c+', 'MarkerSize', 7, 'LineWidth', 1);
        colormap(ax, cmap); xlabel('In-phase'); ylabel('Quadrature');
        title('density, log_{10}(1 + counts)', 'FontSize', 8);
        dark_axes(ax);
        % row 3: phase error
        ax = subplot(3,6,12+ci); hold on;
        if o.snr > 10
            pdf_hist(o.phase_err, linspace(-0.15, 0.15, 81), [0.30 0.75 1.00]); xlim([-0.15 0.15]);
        else
            pdf_hist(o.phase_err, linspace(-pi, pi, 81), [1.00 0.35 0.35]); xlim([-pi pi]);
        end
        xlabel('Phase error [rad]'); ylabel('Density');
        title(sprintf('phase error, \\sigma_\\phi = %.4f rad', o.sigma_phi), 'FontSize', 8);
        dark_axes(ax);
    end

    % column 5: A2 vs pointing state, overlaid with HOLD ON (plain plot, so it can
    % live inside this window; scatterplot() always opens its own figure)
    CF.yrel = [0.01 0.03 0.1 1];               % largest cloud first, smallest on top
    ycol = [1.00 0.35 0.20; 1.00 0.80 0.20; 0.40 0.85 1.00; 0.85 0.85 1.00];
    ax = subplot(3,6,[5 6 11 12]); hold on;
    leg = cell(1, numel(CF.yrel)+1); osym = cell(1, numel(CF.yrel));
    for i = 1:numel(CF.yrel)
        o = rx_metrics(kk_link(P, P.f, P.PT, GQ.A0*CF.yrel(i), CF.mod, MC), MC.Npilot);
        osym{i} = o.sym;
        scatter(real(o.sym), imag(o.sym), 3, ycol(i,:), 'filled');
        leg{i} = sprintf('y/A_0 = %g (SNR %.1f dB)', CF.yrel(i), 10*log10(o.snr));
    end
    scatter(real(pts_n), imag(pts_n), 80, [1 1 1], '+', 'LineWidth', 1.2);
    leg{end} = 'ideal';
    axis equal; xlim(CF.lim*[-1 1]); ylim(CF.lim*[-1 1]);
    xlabel('In-phase'); ylabel('Quadrature');
    title(sprintf('A2 (Tx carrier) vs pointing state, CSPR %.0f dB', P.CSPR_dB), 'FontSize', 9);
    lg = legend(leg, 'Location','southoutside', 'FontSize', 8);
    set(lg, 'TextColor','w', 'Color','k', 'EdgeColor',[0.5 0.5 0.5]);
    dark_axes(ax);

    title_compat(sprintf('%s constellation study at y = A_0: effect of carrier/LO placement (colour = sent symbol)', CF.mod), 'w');
    save_fig(PLOT, 'fig36_constellation_study');

    % ---- optional: the same pointing overlay with the Communications Toolbox
    %      scatterplot() + hold on. scatterplot() always opens its own window,
    %      so it cannot sit inside the grid above; it gets a second black window.
    if PLOT.use_scatterplot && exist('scatterplot', 'file') == 2
        h = [];
        for i = 1:numel(CF.yrel), h = cplot(h, osym{i}, ycol(i,:), PLOT); end
        h = cplot(h, pts_n, [1 1 1], PLOT, '+');
        hold off; axis equal; xlim(CF.lim*[-1 1]); ylim(CF.lim*[-1 1]);
        lg = legend(leg, 'Location','southoutside', 'FontSize', 8);
        set(lg, 'TextColor','w', 'Color','k', 'EdgeColor',[0.5 0.5 0.5]);
        title(sprintf('A2 vs pointing state (scatterplot), CSPR %.0f dB', P.CSPR_dB));
        dark_axes(gca);
        save_fig(PLOT, 'fig36d_scatterplot_overlay');
    end
end

%% ===================== 9. OPTIONAL SWEEPS ==================================
if RUN.cspr
    figure('Name','Rate vs CSPR'); hold on; grid on; cl = lines(numel(MODS));
    for m = 1:numel(MODS)
        Ptmp = P;   % copy so PAPR struct is carried along
        R = arrayfun(@(c) R_erg(setf(Ptmp,'CSPR_dB',c), P.f, P.PT, MODS{m}, 'hard'), SW.CSPR_dB);
        plot(SW.CSPR_dB, R/1e9, 'LineWidth',1.5, 'Color',cl(m,:), 'DisplayName', MODS{m});
        % mark PAPR threshold as vertical dashed line
        xp = 10*log10(P.PAPR.(MODS{m}));
        if exist('xline', 'file')
            xline(xp, '--', MODS{m}, 'Color',cl(m,:), 'HandleVisibility','off', 'FontSize',7);
        else
            plot([xp xp], [0 2*P.B/1e9], '--', 'Color',cl(m,:), 'HandleVisibility','off');
        end
    end
    xlabel('CSPR [dB]'); ylabel('Ergodic rate [Gbps]'); legend('Location','best');
    title({sprintf('Analytic rate vs CSPR (dashed = PAPR, 10^{-4} CCDF), pre-amp %s', onoff(P.pre_on)), ...
           'noise model verified for \kappa \geq 12 dB only; use RUN.cspr\_wave for the optimum'}, 'FontSize', 9);
    save_fig(PLOT, 'fig_opt_rate_vs_cspr');
end
if RUN.jitter
    for m = 1:numel(MODS)
        figure('Name',['Rate vs jitter ' MODS{m}]); hold on; grid on;
        for i = 1:numel(SW.z_vec)
            Pz = setf(P,'z',SW.z_vec(i));
            R = arrayfun(@(s) R_erg(setf(Pz,'sigma_s0',s), P.f, P.PT, MODS{m}, 'hard'), SW.sig0_vec);
            plot(SW.sig0_vec, R/1e9, 'LineWidth',1.5, 'DisplayName', sprintf('z = %d km', SW.z_vec(i)/1e3));
        end
        xlabel('\sigma_{s0} [m]'); ylabel('Ergodic rate [Gbps]'); legend; title(['Rate vs jitter, ' MODS{m}]);
        save_fig(PLOT, ['fig_opt_rate_vs_jitter_' MODS{m}]);
    end
end
if RUN.outage
    figure('Name','Outage rate'); hold on; grid on;
    for m = 1:numel(MODS)
        R = arrayfun(@(e) R_out(P, P.f, P.PT, MODS{m}, e), SW.eps_vec);
        semilogx(SW.eps_vec, R/1e9, 'LineWidth',1.5, 'DisplayName', MODS{m});
    end
    set(gca,'XScale','log'); xlabel('\epsilon'); ylabel('R_\epsilon [Gbps]'); legend; title('\epsilon-outage rate');
    save_fig(PLOT, 'fig_opt_outage_rate');
end
if RUN.cspr_wave
    % Waveform SNR vs CSPR (solid) against the analytic model (dashed; recon term only if P.recon).
    % This is the reference for the optimum CSPR: the waveform includes every KK effect.
    figure('Name','SNR vs CSPR (waveform)','Color','w'); hold on; grid on;
    cl = lines(numel(SW.cspr_yrel)); lg = {};
    for i = 1:numel(SW.cspr_yrel)
        y = GQ.A0*SW.cspr_yrel(i); gw = zeros(size(SW.CSPR_wave_dB)); ga = gw;
        for j = 1:numel(SW.CSPR_wave_dB)
            Pk = setf(P, 'CSPR_dB', SW.CSPR_wave_dB(j));
            gw(j) = kk_wave_mc(Pk, P.f, P.PT, y, MC.wave_mod, MC).snr;
            ga(j) = rx_snr(Pk, P.f, P.hPL*P.PT*y, MC.wave_mod);
        end
        plot(SW.CSPR_wave_dB, 10*log10(gw), '-o', 'Color', cl(i,:), 'LineWidth', 1.5);
        plot(SW.CSPR_wave_dB, 10*log10(ga), '--', 'Color', cl(i,:), 'LineWidth', 1.2);
        lg = [lg, {sprintf('waveform, y/A_0 = %g', SW.cspr_yrel(i)), sprintf('analytic, y/A_0 = %g', SW.cspr_yrel(i))}];
        [gb, jb] = max(gw);
        fprintf('Waveform CSPR sweep %s, y/A0 = %g: best CSPR %d dB (SNR %.1f dB) vs %.1f dB at %d dB\n', ...
            MC.wave_mod, SW.cspr_yrel(i), SW.CSPR_wave_dB(jb), 10*log10(gb), ...
            10*log10(gw(SW.CSPR_wave_dB == P.CSPR_dB)), P.CSPR_dB);
    end
    xlabel('CSPR \kappa [dB]'); ylabel('SNR [dB]'); legend(lg, 'Location','eastoutside');
    title(sprintf('%s: waveform vs analytic SNR across CSPR', MC.wave_mod));
    save_fig(PLOT, 'fig_opt_snr_vs_cspr_waveform');
end
if RUN.wd
    % Fig 3.8: ergodic rate vs detector (aperture) radius
    figure('Name','Rate vs detector radius'); hold on; grid on; cl = lines(numel(MODS));
    fprintf('\nErgodic rate (hard) vs detector radius:\n');
    for m = 1:numel(MODS)
        R = arrayfun(@(w) R_erg(setf(P,'wd',w), P.f, P.PT, MODS{m}, 'hard'), SW.wd_vec);
        plot(SW.wd_vec*100, R/1e9, '-o', 'LineWidth',1.5, 'Color',cl(m,:), 'DisplayName', MODS{m});
        fprintf('  %-4s', MODS{m}); fprintf(' %.2f', R/1e9); fprintf('  Gb/s at w_d = %s cm\n', mat2str(SW.wd_vec*100));
    end
    yl = ylim; plot(P.wd*100*[1 1], yl, 'k--', 'HandleVisibility','off');
    xlabel('Detector radius w_d [cm]'); ylabel('Ergodic rate [Gbps]'); legend('Location','best');
    title(sprintf('Ergodic rate vs detector radius (dashed: default %.0f cm)', P.wd*100));
    save_fig(PLOT, 'fig38_rate_vs_wd');
end
if RUN.itu
    % Fig 3.9: ergodic QPSK rate vs PT at the two ITU-T C-band edges
    figure('Name','Rate vs PT, C-band edges'); hold on; grid on;
    lam = [1530e-9 1565e-9]; st = {'-o','-s'};
    for i = 1:2
        fi = P.c/lam(i); Gi = geom(P, fi, P.z);
        R = arrayfun(@(pt) R_erg(P, fi, pt, 'QPSK', 'hard'), SW.PT_vec);
        plot(SW.PT_vec, R/1e9, st{i}, 'LineWidth',1.5, 'DisplayName', ...
            sprintf('%.0f nm (A_0 = %.3g, \\alpha = %.3f)', lam(i)*1e9, Gi.A0, Gi.alpha));
        fprintf('QPSK at %.0f nm: R(PT = %.2f W) = %.2f Gb/s\n', lam(i)*1e9, P.PT, R_erg(P, fi, P.PT, 'QPSK', 'hard')/1e9);
    end
    xlabel('Transmit power P_T [W]'); ylabel('Ergodic QPSK rate [Gbps]'); legend('Location','southeast');
    title('Ergodic rate at the ITU-T C-band edges');
    save_fig(PLOT, 'fig39_rate_vs_PT_cband');
end
if RUN.csv
    fid = fopen(csv_file,'w');
    fprintf(fid,'Modulation,Freq_THz,PT_W,CSPR_dB,PreAmp,SNRpeak_dB,Rate_hard_Gbps,Rate_soft_Gbps\n');
    for m = 1:numel(MODS)
        for j = 1:numel(SW.f_vec)
            G = geom(P, SW.f_vec(j), P.z);
            for i = 1:numel(SW.PT_vec)
                fprintf(fid,'%s,%.4f,%.4f,%.1f,%d,%.3f,%.6f,%.6f\n', MODS{m}, SW.f_vec(j)/1e12, SW.PT_vec(i), ...
                    P.CSPR_dB, P.pre_on, 10*log10(rx_snr(P,SW.f_vec(j),P.hPL*SW.PT_vec(i)*G.A0)), ...
                    R_erg(P,SW.f_vec(j),SW.PT_vec(i),MODS{m},'hard')/1e9, R_erg(P,SW.f_vec(j),SW.PT_vec(i),MODS{m},'soft')/1e9);
            end
        end
    end
    fclose(fid); fprintf('CSV written: %s\n', csv_file);
end
fprintf('Simulation complete.\n');

%% ===================== 10. CKK / GKK / MP-KK RECEIVER COMPARISON ===========
%  Follows Wen et al., Opt. Express 31(21), 34800 (2023), adapted to the OISL link.
%  An optional front-end distortion (AC-coupled TIA / driver / ADC) acts on the
%  photocurrent before the KK algorithm, built as in the paper's Fig. 3:
%     nonlinearity   v[n] = u[n] + a2*u[n]^2                            (eq. 1)
%     memory         w[n] = v[n] + b*v[n-d] + b*v[n-2d]                 (eq. 2)
%  with a2 and b fixed by the distortion strengths at y = A0
%     S_NL = 10 log10(S/N_NL),   S_ME = 10 log10(S/N_ME)                (eqs. 13-14)
%  Receivers (paper Fig. 7):
%     CKK      : KK directly on the distorted photocurrent
%     CKK + EQ : CKK followed by a symbol-rate FFE after KK (MP.eq_taps taps)
%     GKK      : memoryless polynomial z = d0 + sum_k d_k u^k before KK
%     MP-KK    : memory polynomial z = d0 + sum_k sum_m d_km u[n-m]^k before KK (eq. 11)
%  Compensators are trained on the pilot block against the photocurrent the known
%  pilot symbols predict, by least squares (= converged LMS) or NLMS (Table 1).
%  Metrics: SNR, BER, and the photocurrent NMSE vs the undistorted reference (eq. 15).
%  With MP.SNL_dB = MP.SME_dB = Inf (no distortion) all receivers coincide.
if RUN.mpkk
    MP.mod      = 'QPSK';
    MP.SNL_dB   = 10;            % [dB] nonlinearity strength at y = A0
    MP.SME_dB   = 10;            % [dB] memory-effect strength
    MP.tap_sym  = 0.25;          % [symbols] spacing d of the distortion memory taps
    MP.Kg       = 3;             % GKK polynomial order
    MP.K        = 3;             % MP-KK polynomial order
    MP.M        = 7;             % MP-KK memory length [samples at MC.Nos per symbol]
    MP.eq_taps  = 7;             % CKK + EQ: post-KK FFE taps
    MP.solver   = 'ls';          % 'ls' or 'lms' (normalised LMS, paper Table 1)
    MP.mu       = 0.02;  MP.epochs = 10;   % NLMS step size and passes over the pilots
    MP.yrel     = [1 0.5 0.2 0.1 0.05 0.03 0.02 0.01];   % pointing states for Fig 3.10
    MP.S_sweep  = 4:2:16;        % [dB] S_NL = S_ME sweep (Fig 3.11)
    MP.y_sweep  = 1;             % [y/A0] operating point of the strength and K/M sweeps
    MP.K_sweep  = 1:5;  MP.M_sweep = [1 2 3 5 7 9 13];

    tap  = max(1, round(MP.tap_sym*MC.Nos));
    Lcal = kk_link(P, P.f, P.PT, GQ.A0, MP.mod, MC);     % undistorted on-axis photocurrent
    D    = fe_calibrate(Lcal.Iref, MP.SNL_dB, MP.SME_dB, tap);
    fprintf('\n--- CKK / GKK / MP-KK (%s): S_NL = %.1f dB (a2 = %.3f), S_ME = %.1f dB (b = %.3f, d = %d samples) ---\n', ...
        MP.mod, D.SNL_dB, D.a2, D.SME_dB, D.b, D.tap);
    rxs = mpkk_receivers(D, MP);
    nR  = numel(rxs.opt);

    % ---- Fig 3.10: SNR and BER vs pointing state ----
    snr = zeros(nR, numel(MP.yrel)); ber = snr; nerr = snr; nmse = snr;
    for i = 1:numel(MP.yrel)
        for j = 1:nR
            L = kk_link(P, P.f, P.PT, GQ.A0*MP.yrel(i), MP.mod, MC, rxs.opt{j});
            o = rx_metrics(L, MC.Npilot);
            snr(j,i) = o.snr; ber(j,i) = o.ber; nerr(j,i) = o.nerr; nmse(j,i) = L.nmse;
        end
    end
    fprintf('%-24s', 'y/A0'); fprintf('%8g', MP.yrel); fprintf('\n');
    for j = 1:nR
        fprintf('%-24s', [rxs.name{j} ' SNR[dB]']); fprintf('%8.1f', 10*log10(snr(j,:))); fprintf('\n');
    end
    for j = 1:nR
        fprintf('%-24s', [rxs.name{j} ' BER']); fprintf('%8.1e', ber(j,:)); fprintf('\n');
    end
    figure('Name','CKK GKK MP-KK vs pointing','Color','w','Position',[100 100 1000 420]);
    subplot(1,2,1); hold on; grid on; box on;
    for j = 1:nR
        semilogx(MP.yrel, 10*log10(snr(j,:)), rxs.style{j}, 'Color', rxs.col(j,:), 'LineWidth', 1.4);
    end
    set(gca,'XScale','log'); xlabel('y / A_0'); ylabel('SNR [dB]');
    legend(rxs.name, 'Location','southeast'); title('SNR');
    subplot(1,2,2); hold on; grid on; box on;
    for j = 1:nR      % lines: BER predicted from the measured SNR; markers: counted errors
        loglog(MP.yrel, max(ber_theory(MP.mod, snr(j,:)), 1e-12), rxs.style{j}(1:end-1), ...
            'Color', rxs.col(j,:), 'LineWidth', 1.2);
    end
    for j = 1:nR
        k = nerr(j,:) > 0;
        loglog(MP.yrel(k), ber(j,k), rxs.style{j}(end), 'Color', rxs.col(j,:), ...
            'MarkerSize', 7, 'HandleVisibility','off');
    end
    set(gca,'XScale','log','YScale','log'); ylim([1e-8 0.5]);
    plot(MP.yrel([1 end]), [3.8e-3 3.8e-3], 'k:', 'HandleVisibility','off');
    xlabel('y / A_0'); ylabel('BER');
    title('BER (lines: from SNR; markers: counted; dotted: HD-FEC)');
    title_compat(sprintf('%s receivers with front-end distortion, S_{NL} = %.0f dB, S_{ME} = %.0f dB', ...
        MP.mod, D.SNL_dB, D.SME_dB));
    save_fig(PLOT, 'fig310_mpkk_vs_pointing');

    % ---- Fig 3.11: SNR and photocurrent NMSE vs distortion strength (paper Fig. 8b) ----
    yS = GQ.A0*MP.y_sweep;
    snrS = zeros(nR, numel(MP.S_sweep)); nmseS = snrS;
    for i = 1:numel(MP.S_sweep)
        Di = fe_calibrate(Lcal.Iref, MP.S_sweep(i), MP.S_sweep(i), tap);
        ri = mpkk_receivers(Di, MP);
        for j = 1:nR
            L = kk_link(P, P.f, P.PT, yS, MP.mod, MC, ri.opt{j});
            snrS(j,i) = rx_metrics(L, MC.Npilot).snr; nmseS(j,i) = L.nmse;
        end
    end
    figure('Name','MP-KK vs distortion strength','Color','w','Position',[120 120 1000 420]);
    subplot(1,2,1); hold on; grid on; box on;
    for j = 1:nR, plot(MP.S_sweep, 10*log10(snrS(j,:)), rxs.style{j}, 'Color', rxs.col(j,:), 'LineWidth', 1.4); end
    xlabel('S_{NL} = S_{ME} [dB]'); ylabel('SNR [dB]'); legend(rxs.name, 'Location','southeast');
    title(sprintf('SNR at y/A_0 = %g', MP.y_sweep));
    subplot(1,2,2); hold on; grid on; box on;
    for j = 1:nR, plot(MP.S_sweep, nmseS(j,:), rxs.style{j}, 'Color', rxs.col(j,:), 'LineWidth', 1.4); end
    xlabel('S_{NL} = S_{ME} [dB]'); ylabel('Photocurrent NMSE [dB]');
    title('Error vs undistorted photocurrent (eq. 15)');
    save_fig(PLOT, 'fig311_mpkk_vs_strength');

    % ---- Fig 3.12: MP-KK vs polynomial order K and memory length M (paper Fig. 6) ----
    snrK = zeros(size(MP.K_sweep)); snrM = zeros(size(MP.M_sweep));
    for i = 1:numel(MP.K_sweep)
        c = struct('type','mpkk','K',MP.K_sweep(i),'M',MP.M,'solver',MP.solver,'mu',MP.mu,'epochs',MP.epochs);
        snrK(i) = rx_metrics(kk_link(P, P.f, P.PT, yS, MP.mod, MC, struct('dist',D,'comp',c)), MC.Npilot).snr;
    end
    for i = 1:numel(MP.M_sweep)
        c = struct('type','mpkk','K',MP.K,'M',MP.M_sweep(i),'solver',MP.solver,'mu',MP.mu,'epochs',MP.epochs);
        snrM(i) = rx_metrics(kk_link(P, P.f, P.PT, yS, MP.mod, MC, struct('dist',D,'comp',c)), MC.Npilot).snr;
    end
    fprintf('MP-KK vs K (M = %d):', MP.M); fprintf(' K=%d %.1f dB', [MP.K_sweep; 10*log10(snrK)]); fprintf('\n');
    fprintf('MP-KK vs M (K = %d):', MP.K); fprintf(' M=%d %.1f dB', [MP.M_sweep; 10*log10(snrM)]); fprintf('\n');
    figure('Name','MP-KK order and memory','Color','w','Position',[140 140 1000 380]);
    subplot(1,2,1); plot(MP.K_sweep, 10*log10(snrK), '-o', 'LineWidth', 1.4); grid on;
    xlabel('Polynomial order K'); ylabel('SNR [dB]'); title(sprintf('MP-KK, M = %d', MP.M));
    subplot(1,2,2); plot(MP.M_sweep, 10*log10(snrM), '-s', 'LineWidth', 1.4); grid on;
    xlabel('Memory length M [samples]'); ylabel('SNR [dB]'); title(sprintf('MP-KK, K = %d', MP.K));
    save_fig(PLOT, 'fig312_mpkk_K_M');

    % ---- Fig 3.13: constellation overlay, CKK vs MP-KK, with HOLD ON ----
    oC = rx_metrics(kk_link(P, P.f, P.PT, yS, MP.mod, MC, rxs.opt{2}), MC.Npilot);
    oM = rx_metrics(kk_link(P, P.f, P.PT, yS, MP.mod, MC, rxs.opt{end}), MC.Npilot);
    [pts_c, ~] = const_map(MP.mod);
    h = [];
    h = cplot(h, oC.sym, rxs.col(2,:),   PLOT);
    h = cplot(h, oM.sym, rxs.col(end,:), PLOT);
    h = cplot(h, pts_c/sqrt(mean(abs(pts_c).^2)), [0 0 0], PLOT, 'k+');
    hold off; axis equal; grid on; xlim(1.8*[-1 1]); ylim(1.8*[-1 1]);
    legend({sprintf('CKK (SNR %.1f dB)', 10*log10(oC.snr)), sprintf('MP-KK (SNR %.1f dB)', 10*log10(oM.snr)), 'ideal'}, ...
        'Location','southoutside');
    title(sprintf('%s at y/A_0 = %g, S_{NL} = S_{ME} = %.0f dB', MP.mod, MP.y_sweep, MP.SNL_dB));
    save_fig(PLOT, 'fig313_ckk_vs_mpkk_constellation');
end

%% ===================== 11. DOPPLER RESIDUAL AND RECOVERY ===================
%  With a transmitter carrier the common Doppler shift f*v/c cancels in detection,
%  but the data sideband sits f_sh above the carrier and the symbol clock is scaled,
%  so a residual offset f_sh*v/c (about 263 kHz at 15 km/s) and a v/c (50 ppm) clock
%  scaling remain. This section shows (a) what they do without recovery and (b) that
%  a pilot-based frequency/clock estimate plus blind phase tracking removes them.
%  Frame: MC.Npilot pilots at the start and MC.Npost at the end.
%  Not modelled: Tx/Rx oscillator mismatch (needs timing recovery in any receiver).
if RUN.doppler
    DP.yrel = [1 0.3 0.1 0.03 0.0123 0.01];   % 0.0123*A0 = QPSK HD-FEC threshold
    DP.v    = SW.v_max;
    cases   = {'no Doppler',               0,     false;
               'Doppler, no recovery',     DP.v,  false;
               'Doppler + recovery',       DP.v,  true;
               'Doppler (-v) + recovery', -DP.v,  true};
    snr = zeros(size(cases,1), numel(DP.yrel)); dfh = snr;
    for i = 1:size(cases,1)
        for j = 1:numel(DP.yrel)
            L = kk_link(P, P.f, P.PT, GQ.A0*DP.yrel(j), MC.wave_mod, MC, ...
                        struct('v', cases{i,2}, 'track', cases{i,3}));
            snr(i,j) = rx_metrics(L, MC.Npilot).snr; dfh(i,j) = L.df_hat;
        end
    end
    fprintf('\n--- Doppler residual, %s, v = %.0f km/s (residual offset f_sh*v/c = %.0f kHz, clock %.0f ppm) ---\n', ...
        MC.wave_mod, DP.v/1e3, (1+MC.beta)*P.B/2*DP.v/P.c/1e3, DP.v/P.c*1e6);
    fprintf('%-26s', 'y/A0'); fprintf('%9g', DP.yrel); fprintf('\n');
    for i = 1:size(cases,1)
        fprintf('%-26s', [cases{i,1} ' [dB]']); fprintf('%9.2f', 10*log10(snr(i,:))); fprintf('\n');
    end
    fprintf('%-26s', 'estimated offset [kHz]'); fprintf('%9.1f', dfh(3,:)/1e3); fprintf('\n');
    figure('Name','Doppler residual','Color','w'); hold on; grid on; box on;
    st = {'k--','r-o','b-s','c-^'};
    for i = 1:size(cases,1)
        semilogx(DP.yrel, 10*log10(snr(i,:)), st{i}, 'LineWidth', 1.4);
    end
    set(gca, 'XScale', 'log'); xlabel('y / A_0'); ylabel('SNR [dB]');
    legend(cases(:,1), 'Location', 'southeast');
    title(sprintf('A2 with residual Doppler (v = %.0f km/s): f_{sh}v/c and v/c clock scaling', DP.v/1e3));
    save_fig(PLOT, 'fig_doppler_recovery');
end

%% ===================== 12. RECEIVER-LO COMPARISON (fair B1) ==================
%  A2 (transmitter carrier, with Doppler and recovery) against B1 (receiver LO)
%  with an ideal LO and with a realistic LO: residual frequency offset after
%  locking MC.lo_real_df and combined Tx + LO linewidth MC.lo_real_lw. B1 gets the
%  same frequency and phase recovery, with the phase-tracking window chosen per
%  case from SW.trk_win_vec to give B1 its best result (optimistic for B1).
if RUN.lo_compare
    LC.yrel = [1 0.3 0.1 0.03];
    LC.lw   = MC.lo_real_lw;
    names = [{'A2 Tx carrier (Doppler + recovery)', 'B1 ideal LO'}, ...
             arrayfun(@(w) sprintf('B1 LO: %g MHz offset, %g kHz linewidth', ...
                     MC.lo_real_df/1e6, w/1e3), LC.lw, 'UniformOutput', false)];
    snr = zeros(numel(names), numel(LC.yrel)); wbest = snr;
    for j = 1:numel(LC.yrel)
        y = GQ.A0*LC.yrel(j);
        snr(1,j) = rx_metrics(kk_link(P, P.f, P.PT, y, MC.wave_mod, MC, ...
                   struct('config','tx_carrier','v',SW.v_max,'track',true)), MC.Npilot).snr;
        snr(2,j) = rx_metrics(kk_link(P, P.f, P.PT, y, MC.wave_mod, MC, ...
                   struct('config','rx_lo')), MC.Npilot).snr;
        for i = 1:numel(LC.lw)
            best = 0;
            for W = SW.trk_win_vec
                M2 = MC; M2.lo_df = MC.lo_real_df; M2.lo_lw = LC.lw(i); M2.trk_win = W;
                o = rx_metrics(kk_link(P, P.f, P.PT, y, MC.wave_mod, M2, ...
                    struct('config','rx_lo','track',true)), MC.Npilot);
                if o.snr > best, best = o.snr; wbest(2+i,j) = W; end
            end
            snr(2+i,j) = best;
        end
    end
    fprintf('\n--- Receiver-LO comparison, %s (B1: no EDFA, LO %.2f mW) ---\n', MC.wave_mod, ...
        MC.lo_kappa*P.hPL*MC.lo_yref*P.PT*1e3);
    fprintf('%-44s', 'y/A0'); fprintf('%9g', LC.yrel); fprintf('\n');
    for i = 1:numel(names)
        fprintf('%-44s', [names{i} ' [dB]']); fprintf('%9.2f', 10*log10(snr(i,:))); fprintf('\n');
    end
    for i = 1:numel(LC.lw)
        fprintf('%-44s', sprintf('  best tracking window, %g kHz [symbols]', LC.lw(i)/1e3));
        fprintf('%9d', wbest(2+i,:)); fprintf('\n');
    end
    figure('Name','Receiver-LO comparison','Color','w'); hold on; grid on; box on;
    st = {'k-o','b--','r-s','m-^','g-d'};
    for i = 1:numel(names)
        semilogx(LC.yrel, 10*log10(snr(i,:)), st{i}, 'LineWidth', 1.4);
    end
    set(gca, 'XScale', 'log'); xlabel('y / A_0'); ylabel('SNR [dB]');
    legend(names, 'Location', 'southeast', 'FontSize', 8);
    title('Transmitter carrier (A2) vs receiver LO (B1), same recovery DSP');
    save_fig(PLOT, 'fig_lo_comparison');
end

%% ===================== 13. REALISTIC POINTING: POWER AND RANGE ===============
%  With the stress-case jitter (sigma_s0 = 2 m) the link is below sensitivity 41 %
%  of the time. For realistic jitter the useful question is the transmit power (or
%  range) needed for a target availability. With the power-law CDF of y:
%     P(P_rx < P_sens) = (P_sens / (h_PL P_T A0))^alpha
%  so the minimum transmit power for an outage probability eps is
%     P_T,min = P_sens / (h_PL A0 eps^(1/alpha))
%  (analytic model, kappa = 12 dB). Note: the jitter model of Shang et al. gives
%  sigma_s in metres; for an angular jitter (sigma_s = z*sigma_theta) alpha would
%  not depend on range. See docs/open_questions.md.
if RUN.pointing
    PO.sig0 = SW.sig0_real;
    PO.eps  = SW.avail_eps;
    PO.mod  = 'QPSK';
    Psens   = sensitivity(P, P.f, PO.mod);
    Psens_n = sensitivity(setf(P,'pre_on',false), P.f, PO.mod);
    PTmin = zeros(numel(PO.eps), numel(PO.sig0)); PTmin_n = PTmin;
    al = zeros(size(PO.sig0)); Rerg = al;
    for i = 1:numel(PO.sig0)
        Pi = setf(P, 'sigma_s0', PO.sig0(i)); Gi = geom(Pi, P.f, P.z); al(i) = Gi.alpha;
        PTmin(:,i)   = Psens  ./ (P.hPL*Gi.A0*PO.eps(:).^(1/Gi.alpha));
        PTmin_n(:,i) = Psens_n./ (P.hPL*Gi.A0*PO.eps(:).^(1/Gi.alpha));
        Rerg(i) = R_erg(Pi, P.f, P.PT, PO.mod, 'hard');
    end
    fprintf('\n--- Realistic pointing (%s, z = %.0f km, kappa = %.0f dB) ---\n', PO.mod, P.z/1e3, P.CSPR_dB);
    fprintf('Sensitivity: %.1f dBm (EDFA), %.1f dBm (photodiode only)\n', ...
        10*log10(Psens/1e-3), 10*log10(Psens_n/1e-3));
    fprintf('%-10s %-7s %-14s', 'sigma_s0', 'alpha', 'R_erg(0.3 W)');
    fprintf(' P_T,min @%.0e (EDFA / no EDFA)', PO.eps); fprintf('\n');
    for i = 1:numel(PO.sig0)
        fprintf('%-10.2f %-7.2f %-14s', PO.sig0(i), al(i), sprintf('%.2f Gb/s', Rerg(i)/1e9));
        fprintf('   %7.1f / %7.1f dBm       ', [10*log10(PTmin(:,i)/1e-3) 10*log10(PTmin_n(:,i)/1e-3)].');
        fprintf('\n');
    end
    % feasible range window at P.PT for each availability target. In this jitter
    % model sigma_s is a length that grows slowly with z while the beam radius grows
    % in proportion to z, so outage rises at short AND at long range.
    zg = SW.z_range; Rmin = NaN(numel(PO.eps), numel(PO.sig0)); Rmax = Rmin;
    for i = 1:numel(PO.sig0)
        Pi = setf(P, 'sigma_s0', PO.sig0(i));
        pout = zeros(size(zg));
        for k = 1:numel(zg)
            Gk = geom(Pi, P.f, zg(k));
            pout(k) = min(1, (Psens/(P.hPL*P.PT*Gk.A0))^Gk.alpha);
        end
        for e = 1:numel(PO.eps)
            ok = find(pout <= PO.eps(e));
            if ~isempty(ok), Rmin(e,i) = zg(ok(1)); Rmax(e,i) = zg(ok(end)); end
        end
    end
    fprintf('Feasible range window at P_T = %.2f W (EDFA), km (NaN = none in %.0f-%.0f km):\n', ...
        P.PT, zg(1)/1e3, zg(end)/1e3);
    for e = 1:numel(PO.eps)
        fprintf('  outage %.0e:', PO.eps(e));
        for i = 1:numel(PO.sig0)
            fprintf('  [%g m] %4.0f-%4.0f', PO.sig0(i), Rmin(e,i)/1e3, Rmax(e,i)/1e3);
        end
        fprintf('\n');
    end
    figure('Name','Realistic pointing','Color','w','Position',[100 100 1000 400]);
    subplot(1,2,1); hold on; grid on; box on; cl = lines(numel(PO.eps));
    for e = 1:numel(PO.eps)
        plot(PO.sig0, 10*log10(PTmin(e,:)/1e-3), '-', 'Color', cl(e,:), 'LineWidth', 1.5);
        plot(PO.sig0, 10*log10(PTmin_n(e,:)/1e-3), '--', 'Color', cl(e,:), 'LineWidth', 1.2);
    end
    plot(PO.sig0([1 end]), 10*log10(P.PT/1e-3)*[1 1], 'k:');
    ylim([0 60]); xlabel('\sigma_{s0} [m]'); ylabel('Required P_T [dBm]');
    lg = {}; for e = 1:numel(PO.eps), lg = [lg, {sprintf('outage %.0e, EDFA', PO.eps(e)), ...
                                                   sprintf('outage %.0e, no EDFA', PO.eps(e))}]; end
    legend([lg, {sprintf('P_T = %.1f W', P.PT)}], 'Location', 'northwest', 'FontSize', 7);
    title(sprintf('Transmit power for a target outage, z = %.0f km', P.z/1e3));
    subplot(1,2,2); hold on; grid on; box on;
    for e = 1:numel(PO.eps)
        plot(PO.sig0, Rmax(e,:)/1e3, '-o', 'Color', cl(e,:), 'LineWidth', 1.5);
        plot(PO.sig0, Rmin(e,:)/1e3, '--o', 'Color', cl(e,:), 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    xlabel('\sigma_{s0} [m]'); ylabel('Range [km]');
    legend(arrayfun(@(x) sprintf('outage %.0e (solid max, dashed min)', x), PO.eps, 'UniformOutput', false), ...
        'Location', 'northeast', 'FontSize', 8);
    title(sprintf('Feasible range window at P_T = %.1f W (EDFA)', P.PT));
    save_fig(PLOT, 'fig_pointing_power_range');
end

%% ===================== 14. BER NEAR THE FEC THRESHOLD =======================
%  Counts bit errors over SW.ber_runs x 2^SW.ber_log2N symbols at the pointing
%  states where the analytic BER equals each value in SW.ber_targets, to validate
%  the BER (not only the SNR) in the region that matters for the FEC.
if RUN.berfec
    BF.mod = MC.wave_mod;
    fprintf('\n--- BER near the FEC threshold (%s, %d runs x 2^%d symbols) ---\n', ...
        BF.mod, SW.ber_runs, SW.ber_log2N);
    for bt = SW.ber_targets
        gt = fzero(@(gd) ber_theory(BF.mod, 10^(gd/10)) - bt, [-10 30]);
        Prx = 10^fzero(@(lp) 10*log10(rx_snr(P, P.f, 10^lp, BF.mod)) - gt, [-14 0]);
        y = Prx/(P.hPL*P.PT);
        nerr = 0; nbits = 0; snr_acc = 0;
        for run = 1:SW.ber_runs
            M2 = MC; M2.Nsym = 2^SW.ber_log2N; M2.seed = 100 + run;
            o = rx_metrics(kk_link(P, P.f, P.PT, y, BF.mod, M2), M2.Npilot);
            nerr = nerr + o.nerr; nbits = nbits + o.nbits; snr_acc = snr_acc + o.snr;
        end
        pb = nerr/nbits; ci = 1.96*sqrt(pb*(1-pb)/nbits);
        fprintf('target BER %.1e (g = %.2f dB, y/A0 = %.4f): counted %.3e +- %.1e (%d errors / %d bits), mean SNR %.2f dB\n', ...
            bt, gt, y/GQ.A0, pb, ci, nerr, nbits, 10*log10(snr_acc/SW.ber_runs));
    end
end

%% ===================== 15. FUNCTIONS ======================================
%  Waveform core in separate files (also used by tests/check_oisl_kk.m):
%    kk_link.m, rx_metrics.m, const_map.m, rrc.m, fe_calibrate.m
function G = geom(P, f, z)
    G.lam   = P.c/f;
    G.theta = G.lam/(pi*P.w0);
    G.wz    = z*tan(G.theta);                         % Chapter 2, eq. (2.3) without w0 (far field)
    G.A0    = 2*P.wd^2/G.wz^2;
    G.sig   = P.sigma_s0*exp(P.k0*z/P.d0);
    G.alpha = G.wz^2/(4*G.sig^2);
end

function [g, T] = rx_snr(P, f, Prx, mod)
% Total electrical SNR g (= Es/N0) for received optical power Prx (vector) - eqs. (3.13)-(3.16)
% Optional 4th arg 'mod': adds the heuristic reconstruction penalty T.recon
% (see recon_term) when P.recon is true. Off by default: the waveform simulation
% shows it is far too pessimistic. The validated model is valid for kappa >= 12 dB.
    k    = 10^(P.CSPR_dB/10);  zeta = 1 + 2/k;
    hnu  = P.h*f;
    if P.pre_on
        Gr = 10^(P.G_dB/10); F = 10^(P.NF_dB/10);
        S  = hnu*(F*Gr - 1)/2;          % ASE PSD per polarisation [W/Hz]
    else
        Gr = 1; S = 0;
    end
    Pc  = Gr*k/(1+k)*Prx;              % carrier power at photodiode
    Ps  = Gr/(1+k)*Prx;                % signal power at photodiode
    sig = 2*P.eta^2*Pc.*Ps;            % mean-square cross-beat current [A^2]
    Ibar = P.eta*(Pc + Ps + 2*S*P.Bf);
    RIN = 10^(P.RIN_dB/10);
    T.shot    = zeta*(2*P.q*Ibar)*P.B./sig;
    T.thermal = zeta*(P.i_th^2)*P.B./sig;
    T.RIN     = zeta*(RIN*(P.eta*(Pc+Ps)).^2)*P.B./sig;
    if P.pre_on
        T.ASE    = (1+zeta)*S*P.B./Ps;                              % signal-ASE + carrier-ASE (both sides)
        T.ASExASE = zeta*(4*P.eta^2*S^2*(P.Bf - P.B/2))*P.B./sig;  % ASE-ASE beat, both pols
    end
    % optional heuristic reconstruction penalty
    if nargin >= 4 && ~isempty(mod) && isfield(P,'recon') && P.recon
        T.recon = recon_term(P, mod) + zeros(size(Prx));          % same size as Prx (plots, sums)
    end
    fn = fieldnames(T); tot = 0;
    for i = 1:numel(fn), tot = tot + T.(fn{i}); end
    g = 1./tot;
end

function r = recon_term(P, mod)
% Heuristic KK reconstruction penalty 1/g_rec = PAPR_s * P_viol(kappa), with
% P_viol = exp(-kappa) (QPSK, complex Gaussian) or erfc(sqrt(kappa/2)) (BPSK/OOK).
% NOT validated: 35-50 dB more pessimistic than the waveform simulation.
    k = 10^(P.CSPR_dB/10);
    switch mod
        case 'QPSK',          pviol = exp(-k);
        case {'BPSK','OOK'},  pviol = erfc(sqrt(k/2));
    end
    r = P.PAPR.(mod) * max(pviol, 1e-30);
end

function s = onoff_str(b, a, c), if b, s = a; else, s = c; end, end

function m = mi(mod, g, kind)
% Information bits per symbol. kind = 'hard' (BSC per bit) or 'soft' (GMI)
    qf = @(x) 0.5*erfc(x/sqrt(2));
    h2 = @(p) -p.*log2(max(p,realmin)) - (1-p).*log2(max(1-p,realmin));
    if strcmp(kind,'hard')
        switch mod
            case 'BPSK', m = 1 - h2(qf(sqrt(2*g)));
            case 'OOK',  m = 1 - h2(qf(sqrt(g)));
            case 'QPSK', m = 2*(1 - h2(qf(sqrt(g))));
        end
    else
        switch mod
            case 'BPSK', m = mi_bpsk_soft(g);
            case 'OOK',  m = mi_bpsk_soft(g/2);
            case 'QPSK', m = 2*mi_bpsk_soft(g/2);
        end
    end
end

function I = mi_bpsk_soft(g)
% Soft-decision MI of BPSK at Es/N0 = g, Gauss-Hermite quadrature (no toolbox needed)
    persistent x w
    if isempty(x)
        n = 80; b = sqrt((1:n-1)/2); J = diag(b,1) + diag(b,-1);
        [V, D] = eig(J); x = diag(D); w = sqrt(pi)*V(1,:)'.^2;
    end
    I = zeros(size(g));
    for i = 1:numel(g)
        if g(i) <= 0, I(i) = 0; continue; end
        s2 = 1/(2*g(i)); yv = 1 + sqrt(2*s2)*x; zz = -2*yv/s2;
        l2 = (max(zz,0) + log1p(exp(-abs(zz))))/log(2);      % log2(1+exp(zz)), stable
        I(i) = 1 - sum(w.*l2)/sqrt(pi);
    end
end

function p = ber_theory(mod, g)
    qf = @(x) 0.5*erfc(x/sqrt(2));
    switch mod
        case 'BPSK', p = qf(sqrt(2*g));
        otherwise,   p = qf(sqrt(g));
    end
end

function R = R_erg(P, f, PT, mod, kind)
% Ergodic rate, eq. (3.20): integrate over t = ln y from ln(y_lo) to ln(A0)
    G = geom(P, f, P.z);
    ylo = max(P.pth/(P.hPL*PT), G.A0*1e-30);       % corrected (dimensionless) threshold
    if ylo >= G.A0, R = 0; return; end
    lA0 = log(G.A0);
    fun = @(t) P.B*mi(mod, rx_snr(P, f, P.hPL*PT*exp(t), mod), kind) .* (G.alpha*exp(G.alpha*(t - lA0)));
    R = integral(fun, log(ylo), lA0, 'RelTol',1e-6, 'AbsTol',1e-3);
end

function R = R_out(P, f, PT, mod, eps_)
% epsilon-outage rate, eq. (3.21)
    G = geom(P, f, P.z);
    y = G.A0*eps_^(1/G.alpha);
    if y < P.pth/(P.hPL*PT), R = 0; return; end
    R = P.B*mi(mod, rx_snr(P, f, P.hPL*PT*y, mod), 'hard');
end

function Ps = sensitivity(P, f, mod)
% Received power at which the uncoded BER equals 3.8e-3 (7% HD-FEC threshold)
    gt = fzero(@(gd) ber_theory(mod, 10^(gd/10)) - 3.8e-3, [-10 30]);
    try
        Ps = 10^fzero(@(lp) 10*log10(rx_snr(P, f, 10^lp, mod)) - gt, [-14 0]);
    catch
        Ps = inf;      % reconstruction penalty prevents reaching FEC threshold at this CSPR
    end
end

function Q = setf(P, name, val), Q = P; Q.(name) = val; end

function s = onoff(b), if b, s = 'on'; else, s = 'off'; end, end

function o = kk_wave_mc(P, f, PT, y, mod, MC)
% End-to-end waveform simulation of the Tx-carrier KK link at pointing fraction y.
% Thin wrapper: shared waveform core + pilot-aided receiver metrics.
    o = rx_metrics(kk_link(P, f, PT, y, mod, MC), MC.Npilot);
end

function s = metrics_str(o, mod)
% Two-line metrics summary (cell array) for a panel title
    if o.evm < 1, evs = sprintf('EVM %.1f%%', 100*o.evm); else, evs = 'EVM >100%'; end
    if o.nerr == 0, bs = sprintf('BER 0/%d', o.nbits); else, bs = sprintf('BER %.1e', o.ber); end
    s = {sprintf('SNR %.1f dB | %s', 10*log10(o.snr), evs), ...
         sprintf('\\sigma_\\phi %.3f rad | %s', o.sigma_phi, bs)};
end

function h = cplot(h, x, col, PLOT, style)
% Overlay complex samples on one black axes with HOLD ON. Uses scatterplot() from
% the Communications Toolbox when PLOT.use_scatterplot is true and it is installed,
% otherwise scatter(). col picks the nearest scatterplot colour letter.
    if nargin < 5, style = '.'; end
    usecomm = PLOT.use_scatterplot && exist('scatterplot', 'file') == 2;
    if usecomm
        cmap = 'rycwgbm';
        [~, ci] = min(sum(abs(col - [1 0 0; 1 1 0; 0 1 1; 1 1 1; 0 1 0; 0 0 1; 1 0 1]), 2));
        st = [cmap(ci) style];
        if isempty(h), h = scatterplot(x, 1, 0, st); else, scatterplot(x, 1, 0, st, h); end
        hold on;
    else
        if isempty(h), h = figure('Name','Constellation overlay'); end
        hold on;
        if strcmp(style, '.')
            scatter(real(x), imag(x), 4, col, 'filled');
        else
            scatter(real(x), imag(x), 70, col, style, 'LineWidth', 1);
        end
    end
    set(h, 'Color', 'k');
    if isprop(h, 'InvertHardcopy'), set(h, 'InvertHardcopy', 'off'); end
    dark_axes(gca);
end

function H = hist2d_counts(x, y, edges)
    n  = numel(edges) - 1;
    ix = floor((x - edges(1))/(edges(end)-edges(1))*n) + 1;
    iy = floor((y - edges(1))/(edges(end)-edges(1))*n) + 1;
    ok = ix >= 1 & ix <= n & iy >= 1 & iy <= n;
    H  = accumarray([iy(ok) ix(ok)], 1, [n n]);
end

function pdf_hist(x, edges, col)
    if exist('histcounts', 'file') || exist('histcounts', 'builtin')
        c = histcounts(x(:), edges);
    else
        c = histc(x(:), edges); c = c(1:end-1).';
    end
    w = edges(2) - edges(1);
    bar(edges(1:end-1) + w/2, c/(sum(c)*w), 1, 'FaceColor', col, 'EdgeColor', 'none');
end

function cm = parula_compat()
    if exist('parula', 'file') || exist('parula', 'builtin'), cm = parula(256); else, cm = jet(256); end
end

function title_compat(txt, col)
    if nargin < 2, col = 'k'; end
    if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
        sgtitle(txt, 'FontSize', 10, 'Color', col);
    else
        annotation('textbox', [0 0.955 1 0.04], 'String', txt, 'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', 'Color', col);
    end
end

function dark_axes(ax)
% Black axes with light text and grid, for the black-background figure
    set(ax, 'Color','k', 'XColor',[0.85 0.85 0.85], 'YColor',[0.85 0.85 0.85], ...
        'GridColor',[0.6 0.6 0.6], 'GridAlpha',0.35, 'Box','on');
    grid(ax, 'on');
    set(get(ax,'Title'), 'Color','w');
    set(get(ax,'XLabel'), 'Color',[0.85 0.85 0.85]);
    set(get(ax,'YLabel'), 'Color',[0.85 0.85 0.85]);
end

function save_fig(PLOT, name)
    if ~PLOT.save, return; end
    if ~exist(PLOT.dir, 'dir'), mkdir(PLOT.dir); end
    if exist('exportgraphics', 'file')            % MATLAB R2020a+: tight, vector PDF
        exportgraphics(gcf, fullfile(PLOT.dir, [name '.png']), 'Resolution', 150);
        exportgraphics(gcf, fullfile(PLOT.dir, [name '.pdf']), 'ContentType', 'vector');
    else
        print(gcf, fullfile(PLOT.dir, [name '.png']), '-dpng', '-r150');
        print(gcf, fullfile(PLOT.dir, [name '.pdf']), '-dpdf');
    end
end

function R = mpkk_receivers(D, MP)
% Receiver set compared in Section 10 (first entry: undistorted reference)
    cg = struct('type','gkk', 'K',MP.Kg, 'M',1,    'solver',MP.solver,'mu',MP.mu,'epochs',MP.epochs);
    cm = struct('type','mpkk','K',MP.K,  'M',MP.M, 'solver',MP.solver,'mu',MP.mu,'epochs',MP.epochs);
    R.opt   = {struct(), struct('dist',D), struct('dist',D,'post_eq',MP.eq_taps), ...
               struct('dist',D,'comp',cg), struct('dist',D,'comp',cm)};
    R.name  = {'No distortion','CKK','CKK + post-KK EQ','GKK','MP-KK'};
    R.style = {'--x','-o','-s','-^','-d'};
    R.col   = [0.5 0.5 0.5; 0.25 0.25 0.25; 0 0.45 0.74; 0.47 0.67 0.19; 0.85 0.20 0.20];
end

function papr = signal_papr(mod, B, beta)
% Peak-to-average power ratio of the RRC-shaped signal, defined as the 10^-4 point
% of the instantaneous-power CCDF (2^17 symbols, 8 samples/symbol, fixed seed).
% A CCDF point is stable; the single largest peak of one waveform is not (it
% depends on the random generator). The caller's RNG state is restored.
    Nsym = 2^17; Nos = 8; N = Nsym*Nos; fs = Nos*B;
    fr = ((0:N-1).' - N*((0:N-1).' >= N/2))*fs/N;
    pts = const_map(mod);
    rs = rng; rng(42);
    up = zeros(N,1); up(1:Nos:end) = pts(randi(numel(pts), Nsym, 1));
    rng(rs);
    s  = ifft(fft(up).*rrc(fr/B, beta));
    p  = sort(abs(s).^2/mean(abs(s).^2));
    papr = p(ceil((1 - 1e-4)*N));                % 10^-4 CCDF point (linear)
end
