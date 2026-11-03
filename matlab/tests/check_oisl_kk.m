% CHECK_OISL_KK  Regression checks for oisl_kk_channel.m (MATLAB or Octave).
%
%   Run from the matlab/ folder:   run tests/check_oisl_kk.m      (about 30 s)
%   Part 1 checks the analytic model (oisl_kk_channel.m) against reference values
%   in docs/results.md. Part 2 checks the waveform simulation (kk_link.m,
%   rx_metrics.m) against the analytic model and against physical expectations.
%   Waveform tolerances allow for MATLAB and Octave using different random
%   generators. A failure means the model changed: update the reference on
%   purpose (and the report numbers with it) or fix the change.

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));   % find oisl_kk_channel.m
nfail = 0;
chk = @(name, ok) fprintf('  [%s] %s\n', char(ok*'PASS' + ~ok*'FAIL'), name);

fprintf('Part 1: analytic model (oisl_kk_channel)\n');

r = oisl_kk_channel();                                   % defaults, QPSK, pre-amp on
ok = abs(r.snr_dB - 31.05) < 0.05;            chk('on-axis SNR = 31.05 dB', ok);            nfail = nfail + ~ok;
ok = abs(r.Prx_dBm - (-9.52)) < 0.01;         chk('received power = -9.52 dBm', ok);        nfail = nfail + ~ok;
ok = abs(r.zeta - 1.126) < 1e-3;              chk('zeta = 1.126', ok);                      nfail = nfail + ~ok;
ok = abs(r.rate_hard/1e9 - 13.18) < 0.02;     chk('QPSK ergodic rate (hard) = 13.18 Gb/s', ok); nfail = nfail + ~ok;
ok = abs(r.rate_soft/1e9 - 13.44) < 0.02;     chk('QPSK ergodic rate (soft) = 13.44 Gb/s', ok); nfail = nfail + ~ok;
ok = abs(r.sens_dBm - (-28.6)) < 0.1;         chk('QPSK sensitivity = -28.6 dBm', ok);      nfail = nfail + ~ok;
ok = abs(r.geom.alpha - 0.2006) < 1e-3;       chk('alpha = 0.2006', ok);                    nfail = nfail + ~ok;

r = oisl_kk_channel('pre_on', false);
ok = abs(r.rate_hard/1e9 - 6.70) < 0.02;      chk('QPSK, no pre-amp = 6.70 Gb/s', ok);      nfail = nfail + ~ok;

r = oisl_kk_channel('mod', 'BPSK');
ok = abs(r.rate_hard/1e9 - 6.84) < 0.02;      chk('BPSK ergodic rate (hard) = 6.84 Gb/s', ok); nfail = nfail + ~ok;
ok = abs(r.snr_dB - 31.05) < 0.05;            chk('BPSK on-axis SNR = 31.05 dB (no heuristic penalty)', ok); nfail = nfail + ~ok;
ok = abs(r.PAPR_dB - 7.45) < 0.25;            chk('BPSK PAPR (1e-4 CCDF) = 7.45 +- 0.25 dB', ok); nfail = nfail + ~ok;
r = oisl_kk_channel();
ok = abs(r.PAPR_dB - 6.35) < 0.25;            chk('QPSK PAPR (1e-4 CCDF) = 6.35 +- 0.25 dB', ok); nfail = nfail + ~ok;
r = oisl_kk_channel('mod', 'BPSK', 'recon', true);
ok = r.snr_dB < 30.5;                         chk('recon = true adds the heuristic penalty', ok); nfail = nfail + ~ok;

r = oisl_kk_channel();
ok = abs(r.p_below_sens - 0.414) < 0.005;     chk('time below QPSK sensitivity = 41.4 %', ok); nfail = nfail + ~ok;
r = oisl_kk_channel('sigma_s0', 0.5);
ok = abs(10*log10(r.PT_min_1e2/1e-3) - 11.9) < 0.1;  chk('sigma_s0 = 0.5 m: P_T for 1e-2 outage = 11.9 dBm', ok); nfail = nfail + ~ok;

r = oisl_kk_channel('wd', 0.10);
ok = abs(r.rate_hard/1e9 - 14.09) < 0.02;     chk('wd = 0.10 m -> 14.09 Gb/s', ok);         nfail = nfail + ~ok;

r = oisl_kk_channel('f', 299792458/1565e-9);
ok = abs(r.rate_hard/1e9 - 13.45) < 0.02;     chk('1565 nm -> 13.45 Gb/s', ok);             nfail = nfail + ~ok;

% the caller's random-number state must not be changed by the function
rng(7); a = rand; rng(7); oisl_kk_channel(); b = rand;
ok = (a == b);                                chk('caller RNG state preserved', ok);        nfail = nfail + ~ok;

% invalid input must raise an error
bad = {{'Pt', 0.2}, {'config', 'rx_lo'}, {'mod', '16QAM'}, {'PT'}};
for i = 1:numel(bad)
    try
        oisl_kk_channel(bad{i}{:}); ok = false;
    catch
        ok = true;
    end
    chk(sprintf('rejects invalid input #%d', i), ok); nfail = nfail + ~ok;
end

% ---------------- Part 2: waveform simulation ----------------
fprintf('Part 2: waveform simulation (kk_link, rx_metrics)\n');
r0 = oisl_kk_channel();  P = r0.params;  A0 = r0.geom.A0;
MC = struct('Nsym', 2^13, 'Nos', 8, 'beta', P.beta, 'seed', 1, 'Npilot', 512, 'guard', 16, ...
            'lo_kappa', 10^(12/10), 'lo_yref', A0, 'lo_df', 0, 'lo_lw', 0, ...
            'trk_win', 129, 'Npost', 128);
for yr = [1 0.1 0.03]
    y  = A0*yr;
    o  = rx_metrics(kk_link(P, P.f, P.PT, y, 'QPSK', MC), MC.Npilot);
    ra = oisl_kk_channel('y_override', y);
    d  = 10*log10(o.snr) - ra.snr_dB;
    ok = abs(d) < 0.3;
    chk(sprintf('QPSK waveform SNR at y/A0 = %g within 0.3 dB of analytic (%+.2f dB)', yr, d), ok); nfail = nfail + ~ok;
end
nz = struct('ase', false, 'shot', false, 'th', false, 'rin', false);
o  = rx_metrics(kk_link(P, P.f, P.PT, A0, 'QPSK', MC, struct('noise', nz)), MC.Npilot);
ok = 10*log10(o.snr) > 80;                    chk(sprintf('noiseless KK reconstruction > 80 dB (%.1f dB)', 10*log10(o.snr)), ok); nfail = nfail + ~ok;
oA = rx_metrics(kk_link(P, P.f, P.PT, A0, 'QPSK', MC, struct('config', 'tx_carrier')), MC.Npilot);
oB = rx_metrics(kk_link(P, P.f, P.PT, A0, 'QPSK', MC, struct('config', 'rx_lo')), MC.Npilot);
ok = oB.snr > oA.snr;                         chk('ideal Rx LO (B1) beats Tx carrier (A2)', ok); nfail = nfail + ~ok;
for c = {'no_carrier', 'post_det_lo'}
    o  = rx_metrics(kk_link(P, P.f, P.PT, A0, 'QPSK', MC, struct('config', c{1})), MC.Npilot);
    ok = o.ber > 0.4;                         chk(sprintf('%s fails (BER %.2f)', c{1}, o.ber), ok); nfail = nfail + ~ok;
end
cm = struct('type', 'mpkk', 'K', 3, 'M', 7, 'solver', 'ls');
o  = rx_metrics(kk_link(P, P.f, P.PT, A0, 'QPSK', MC, struct('comp', cm)), MC.Npilot);
ok = abs(10*log10(o.snr) - 10*log10(oA.snr)) < 0.2;
chk('MP-KK on an undistorted link changes nothing', ok); nfail = nfail + ~ok;

% residual Doppler (v = 15 km/s) with and without recovery
for yr = [1 0.1]
    y  = A0*yr;
    o0 = rx_metrics(kk_link(P, P.f, P.PT, y, 'QPSK', MC), MC.Npilot);
    oN = rx_metrics(kk_link(P, P.f, P.PT, y, 'QPSK', MC, struct('v', 15e3)), MC.Npilot);
    oR = rx_metrics(kk_link(P, P.f, P.PT, y, 'QPSK', MC, struct('v', 15e3, 'track', true)), MC.Npilot);
    ok = 10*log10(oN.snr) < 10;
    chk(sprintf('y/A0 = %g: Doppler without recovery breaks the link (%.1f dB)', yr, 10*log10(oN.snr)), ok); nfail = nfail + ~ok;
    d  = 10*log10(oR.snr) - 10*log10(o0.snr); ok = abs(d) < 0.5;
    chk(sprintf('y/A0 = %g: Doppler + recovery within 0.5 dB of no Doppler (%+.2f dB)', yr, d), ok); nfail = nfail + ~ok;
end
% realistic receiver LO (100 MHz residual offset, 200 kHz linewidth) with recovery
M2 = MC; M2.lo_df = 100e6; M2.lo_lw = 200e3; M2.trk_win = 5;
oB = rx_metrics(kk_link(P, P.f, P.PT, A0, 'QPSK', M2, struct('config', 'rx_lo', 'track', true)), M2.Npilot);
ok = 10*log10(oB.snr) > 33 && 10*log10(oB.snr) < 37.3;
chk(sprintf('realistic Rx LO between A2 and ideal LO (%.1f dB)', 10*log10(oB.snr)), ok); nfail = nfail + ~ok;

if nfail == 0
    fprintf('All checks passed.\n');
else
    error('%d check(s) failed.', nfail);
end

