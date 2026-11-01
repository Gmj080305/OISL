% CHECK_OISL_KK  Regression checks for oisl_kk_channel.m (MATLAB or Octave).
%
%   Run from the matlab/ folder:   run tests/check_oisl_kk.m
%   Each check compares against a reference value from the October 2026 runs
%   (docs/results.md). A failure means the model changed: either update the
%   reference on purpose (and the report numbers with it) or fix the change.
%
%   The waveform simulation in oisl_kk_txcarrier.m is checked by its own console
%   output; the expected values are listed in docs/results.md.

addpath(fullfile(fileparts(mfilename('fullpath')), '..'));   % find oisl_kk_channel.m
nfail = 0;
chk = @(name, ok) fprintf('  [%s] %s\n', char(ok*'PASS' + ~ok*'FAIL'), name);

fprintf('oisl_kk_channel regression checks\n');

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

if nfail == 0
    fprintf('All checks passed.\n');
else
    error('%d check(s) failed.', nfail);
end

