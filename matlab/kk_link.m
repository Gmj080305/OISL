function L = kk_link(P, f, PT, y, mod, MC, opt)
% Shared end-to-end waveform link (Tx -> channel -> [EDFA] -> photodiode -> KK DSP
% -> matched filter). Returns symbol-rate samples r and the sent symbols.
%   Part of the KK-OISL simulation (see oisl_kk_txcarrier.m). Works in MATLAB and GNU Octave.
%   opt.config : 'tx_carrier'  A2  carrier + SSB data share PT, both faded by y
%                'rx_lo'       B1  all of PT carries data; LO added just before the
%                                  diode, no EDFA (the LO provides the gain)
%                'no_carrier'  A1/B3 data only, no carrier anywhere
%                'post_det_lo' B2  data only at the diode; the carrier level is
%                                  added digitally after detection, before KK
%   opt.noise  : struct of logicals .ase .shot .th .rin (all true by default),
%                lets each noise term be switched off to test the analytic model
%   opt.dist   : front-end distortion struct from fe_calibrate (empty = none)
%   opt.comp   : pre-KK compensation: .type 'none' | 'gkk' | 'mpkk', .K, .M,
%                .solver 'ls' | 'lms', .mu, .epochs  (trained on the pilot block)
%   opt.post_eq: number of taps of a post-KK symbol-rate FFE (0 = none)
%   opt.v      : relative radial velocity [m/s] (default 0). Doppler time-scales the
%                data envelope by (1 + v/c) and the SSB offset to f_sh*(1 + v/c).
%                The common carrier shift f*v/c cancels in detection (A2) or is
%                assumed removed by LO locking to within MC.lo_df (B1).
%   opt.track  : true = receiver frequency/clock/phase recovery (see recover()):
%                pilot-based frequency estimate, clock correction from the same
%                estimate (A2 only: carrier offset and clock share one Doppler
%                factor), blind M-th power phase tracking. Default false.
    if nargin < 7, opt = struct(); end
    if ~isfield(opt,'config'), opt.config = 'tx_carrier'; end
    if ~isfield(opt,'noise'),  opt.noise = struct('ase',true,'shot',true,'th',true,'rin',true); end
    if ~isfield(opt,'v'),      opt.v = 0; end
    if ~isfield(opt,'track'),  opt.track = false; end
    rng(MC.seed);
    Ns = MC.Nsym; Nos = MC.Nos; N = Ns*Nos; fs = Nos*P.B;
    t  = (0:N-1).'/fs;
    fr = ((0:N-1).' - N*((0:N-1).' >= N/2))*fs/N;          % FFT frequency grid
    [pts, lab] = const_map(mod);
    idx = randi(numel(pts), Ns, 1); a = pts(idx);
    Hr = rrc(fr/P.B, MC.beta);                             % RRC pulse shaping (cyclic)
    up = zeros(N,1); up(1:Nos:end) = a;
    s  = ifft(fft(up).*Hr); s = s/sqrt(mean(abs(s).^2));   % unit average power
    Bo = (1+MC.beta)*P.B; df = fs/N; fsh = ceil((Bo/2)/df)*df + df;
    k  = 10^(P.CSPR_dB/10);
    dop = opt.v/P.c;                                       % Doppler factor v/c
    if dop ~= 0, s = cyclic_scale(s, 1 + dop); end         % envelope time-scaled
    ssb = s.*exp(1i*2*pi*fsh*(1 + dop)*t);                 % SSB offset scaled too
    pre_on = P.pre_on; ElO = 0; P_LO = 0;
    % ---- optical field at the receiver aperture (after the channel)
    switch opt.config
        case 'tx_carrier'
            E = sqrt(P.hPL*y)*(sqrt(k/(1+k)*PT) + sqrt(PT/(1+k))*ssb);
        case 'rx_lo'
            E = sqrt(P.hPL*y*PT)*ssb;
            pre_on = false;
            P_LO = MC.lo_kappa*P.hPL*MC.lo_yref*PT;        % CSPR = lo_kappa at y = lo_yref
            ph = 2*pi*MC.lo_df*t;                          % residual LO frequency offset [Hz]
            if MC.lo_lw > 0                                % combined Tx + LO linewidth [Hz]
                ph = ph + cumsum(sqrt(2*pi*MC.lo_lw/fs)*randn(N,1));
            end
            ElO = sqrt(P_LO)*exp(-1i*ph);
        case {'no_carrier','post_det_lo'}
            E = sqrt(P.hPL*y*PT)*ssb;
        otherwise, error('Unknown config: %s', opt.config);
    end
    % ---- optional pre-amplifier with ASE (both polarisations) and optical BPF of width Bf
    if pre_on
        Gr = 10^(P.G_dB/10); F = 10^(P.NF_dB/10); S = P.h*f*(F*Gr-1)/2;
        bpf = abs(fr - Bo/2) <= P.Bf/2;
        n1 = ifft(fft(sqrt(S*fs/2)*(randn(N,1)+1i*randn(N,1))).*bpf);
        n2 = ifft(fft(sqrt(S*fs/2)*(randn(N,1)+1i*randn(N,1))).*bpf);
        if ~opt.noise.ase, n1 = 0*n1; n2 = 0*n2; end
        Ex = sqrt(Gr)*E + n1; Ey = n2;
    else
        Gr = 1; Ex = E; Ey = 0;
    end
    Ex = Ex + ElO;                                         % B1 only (ElO = 0 otherwise)
    % ---- photodiode: square law + shot, thermal and RIN (white, one-sided PSD Ni)
    I    = P.eta*(abs(Ex).^2 + abs(Ey).^2);
    Idc  = mean(I);
    Ni   = opt.noise.shot*2*P.q*Idc + opt.noise.th*P.i_th^2 + opt.noise.rin*10^(P.RIN_dB/10)*Idc^2;
    I    = I + sqrt(Ni*fs/2)*randn(N,1);
    I    = real(ifft(fft(I).*(abs(fr) <= Bo)));             % electrical bandwidth Bo
    if strcmp(opt.config,'post_det_lo')                    % B2: carrier term added in DSP
        I = I + P.eta*Gr*P.hPL*y*PT*k;
    end
    % ---- noiseless, undistorted reference photocurrent (what training symbols predict)
    Iref = P.eta*abs(sqrt(Gr)*E + ElO).^2;
    Iref = real(ifft(fft(Iref).*(abs(fr) <= Bo)));
    % ---- optional receiver front-end distortion (Section 10, paper eqs. 1-2)
    if isfield(opt,'dist') && ~isempty(opt.dist), I = fe_distort(I, opt.dist); end
    trn = (MC.guard*Nos+1):((MC.Npilot-MC.guard)*Nos);    % training samples (pilot block)
    pay = (MC.Npilot*Nos+1):N;                              % payload samples
    L.nmse_in = nmse_db(I, Iref, pay);
    % ---- optional pre-KK compensation: GKK (memoryless) or MP-KK (memory polynomial)
    if isfield(opt,'comp') && ~strcmp(opt.comp.type,'none')
        I = fe_compensate(I, Iref, trn, opt.comp);
    end
    L.nmse = nmse_db(I, Iref, pay);                        % photocurrent error vs reference
    % ---- KK DSP (eqs. 3.6-3.11)
    In  = I/P.eta; In = max(In, 1e-9*mean(In));
    phi = real(ifft(-1i*sign(fr).*fft(0.5*log(In))));
    Eh  = sqrt(In).*exp(1i*phi);
    xs  = (Eh - mean(Eh)).*exp(-1i*2*pi*fsh*t);            % remove carrier, then down-shift
    r   = ifft(fft(xs).*Hr); r = r(1:Nos:end);             % matched filter, 1 sample/symbol
    L.df_hat = 0; L.eps_hat = 0;
    if opt.track
        [r, L.df_hat, L.eps_hat] = recover(xs, a, pts, Hr, t, fsh, Nos, P.B, MC, ...
                                           strcmp(opt.config, 'tx_carrier'));
    end
    % ---- optional post-KK symbol-rate equaliser (the paper's "CKK + equaliser")
    if isfield(opt,'post_eq') && opt.post_eq > 0
        r = sym_ffe(r, a, MC.Npilot, opt.post_eq);
    end
    L.a = a; L.idx = idx; L.r = r; L.pts = pts; L.lab = lab;
    L.npost = 0; if opt.track, L.npost = MC.Npost; end      % end pilots, excluded from metrics
    L.I = I; L.Iref = Iref;
    L.Pd_dBm = 10*log10(mean(abs(Ex).^2 + abs(Ey).^2)/1e-3); % optical power on the diode
    L.P_LO = P_LO;
end

%% ---- local helpers: front-end distortion and pre-KK compensation ----
function I = fe_distort(I, D)
% AC-coupled front-end (TIA / driver / ADC) model, paper Fig. 3 in series:
%   nonlinearity  v[n] = u[n] + a2*u[n]^2                       (paper eq. 1, K = 2)
%   memory        w[n] = v[n] + b*v[n-d] + b*v[n-2d]            (paper eq. 2, M = 3)
% u is the AC part of the photocurrent in units of its on-axis rms (D.sig_ref),
% so the nonlinearity weakens in a fade, as it would for a real device.
% The DC level is restored ideally afterwards (KK needs it).
    Idc = mean(I); u = (I - Idc)/D.sig_ref;
    v = u + D.a2*u.^2;
    w = v + D.b*(circshift(v, D.tap) + circshift(v, 2*D.tap));
    I = Idc + D.sig_ref*(w - mean(w));
end

function Io = fe_compensate(I, Iref, trn, C)
% Pre-KK inverse model trained on the pilot block (paper eq. 11):
%   GKK  : z[n] = d0 + sum_{k=1..K} d_k * u[n]^k                       (memoryless)
%   MP-KK: z[n] = d0 + sum_{k=1..K} sum_{m=0..M-1} d_km * u[n-m]^k     (memory polynomial)
% For a real, non-negative intensity the paper's x|x|^k terms are plain powers.
% Coefficients by least squares (= converged LMS) or by normalised LMS (paper Table 1).
    mu0 = mean(I(trn)); s0 = std(I(trn)); u = (I - mu0)/s0;
    if strcmp(C.type, 'gkk'), M = 1; else, M = C.M; end
    X  = mp_basis(u, C.K, M);
    dm = mean(Iref(trn)); ds = std(Iref(trn)); d = (Iref - dm)/ds;
    if isfield(C,'solver') && strcmp(C.solver, 'lms')
        w = zeros(size(X,2), 1);
        for ep = 1:C.epochs
            for n = trn
                x = X(n,:).'; e = d(n) - w.'*x;
                w = w + C.mu*e*x/(1e-9 + x.'*x);
            end
        end
    else
        w = X(trn,:)\d(trn);
    end
    Io = dm + ds*(X*w);
end

function X = mp_basis(u, K, M)
    X = ones(numel(u), 1 + K*M); c = 1;
    for m = 0:M-1
        um = circshift(u, m);
        for k = 1:K, c = c + 1; X(:,c) = um.^k; end
    end
end

function v = nmse_db(I, Iref, idx)
% Normalised mean-square error of the photocurrent vs the undistorted reference
% (paper eq. 15), after removing gain and offset
    A = [ones(numel(idx),1) I(idx)]; c = A\Iref(idx);
    v = 10*log10(mean((Iref(idx) - A*c).^2)/var(Iref(idx)));
end

function r = sym_ffe(r, a, Np, T)
% Post-KK complex symbol-rate FFE with T taps, least squares on the pilots
    c = (T-1)/2; R = zeros(numel(r), T);
    for j = 1:T, R(:,j) = circshift(r, j-1-c); end
    w = R(1:Np,:)\a(1:Np);
    r = R*w;
end

%% ---- Doppler model and receiver recovery ----
function y = cyclic_scale(x, a)
% y(t) = x(a*t) for a cyclic, oversampled signal (spline interpolation)
    N  = numel(x); n = (0:N-1).';
    tq = mod(n*a, N);
    xe = [x; x(1:4)]; ne = (0:N+3).';
    y  = interp1(ne, real(xe), tq, 'spline') + 1i*interp1(ne, imag(xe), tq, 'spline');
end

function [r, df_hat, eps_hat] = recover(xs, a, pts, Hr, t, fsh, Nos, B, MC, clock_from_carrier)
% Receiver frequency, clock and phase recovery on the down-shifted KK field xs.
% Frame: MC.Npilot known symbols at the start, MC.Npost known symbols at the end.
%  1. Coarse frequency: slope of the unwrapped phase of r(k)*conj(a(k)) over the
%     start pilots (least squares).
%  2. Fine frequency: after the coarse frequency and clock correction, phase
%     difference between the start and end pilot blocks, divided by their
%     distance in time. MC.guard symbols at each frame edge are not used.
%  3. Transmitter carrier (A2): the residual offset is f_sh*v/c and the clock is
%     scaled by the same v/c, so eps = df/f_sh corrects both. With a receiver LO
%     (B1) the two are unrelated and only the frequency is corrected.
%  4. Remaining phase: blind M-th power (Viterbi-Viterbi) estimator over a moving
%     window of MC.trk_win symbols, excluding the symbol being corrected so its
%     own noise is not removed (which would inflate the measured SNR); the pilot
%     gain in rx_metrics then removes the M-fold ambiguity.
    T   = 1/B;
    mf  = @(x) subsample(ifft(fft(x).*Hr), Nos);
    Ns  = numel(a); g = MC.guard;                           % keep clear of the frame edges
    p1  = (g+1:MC.Npilot).';
    % coarse frequency from the start pilots
    r0  = mf(xs);
    ph  = unwrap(angle(r0(p1).*conj(a(p1))));
    c   = [ones(numel(p1),1) (p1-1)]\ph;
    df_hat = c(2)/(2*pi*T);
    % fine frequency from start and end pilots, after the coarse clock correction
    if MC.Npost > 0
        p2 = (Ns-MC.Npost+1:Ns-g).';
        r1 = mf(correct(xs, df_hat, fsh, t, clock_from_carrier));
        dphi = angle(sum(r1(p2).*conj(a(p2))) * conj(sum(r1(p1).*conj(a(p1)))));
        df_hat = df_hat + dphi/(2*pi*(mean(p2) - mean(p1))*T);
    end
    [xs, eps_hat] = correct(xs, df_hat, fsh, t, clock_from_carrier);
    r   = mf(xs);
    M   = numel(pts);
    z   = (r/mean(abs(r))).^M;
    W   = MC.trk_win; h = ones(W,1);
    zs  = ifft(fft(z).*fft([h; zeros(numel(z)-W,1)]));     % cyclic moving sum
    zs  = circshift(zs, -floor(W/2));                       % centre the window
    zs  = zs - z;                                           % leave the symbol itself out
    th  = unwrap(angle(zs))/M;                              % (else its own noise is removed)
    r   = r.*exp(-1i*th);
end

function [y, eps_hat] = correct(x, df_hat, fsh, t, clock_from_carrier)
% Remove the frequency offset and, for the transmitter carrier, the clock scaling
    y = x.*exp(-1i*2*pi*df_hat*t);
    eps_hat = 0;
    if clock_from_carrier
        eps_hat = df_hat/fsh;
        y = cyclic_scale(y, 1/(1 + eps_hat));
    end
end

function y = subsample(x, Nos)
    y = x(1:Nos:end);
end
