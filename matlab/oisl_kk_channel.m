function result = oisl_kk_channel(varargin)
% OISL_KK_CHANNEL  Analytical KK self-coherent OISL channel model.
%
%   result = oisl_kk_channel()
%       Run with all default parameters.
%
%   result = oisl_kk_channel('param1', val1, 'param2', val2, ...)
%       Override specific parameters as name-value pairs.
%
%   result = oisl_kk_channel(S)
%       Override parameters from struct S (fields override defaults).
%
%   Computes analytical SNR, noise budget, ergodic/outage data rates, and
%   receiver sensitivity for an optical inter-satellite link (OISL) using a
%   Kramers-Kronig (KK) self-coherent receiver with transmitter-side carrier.
%
%   The analytic model describes the transmitter-side-carrier link (config
%   'tx_carrier') only. The other carrier/LO placements (B1 Rx LO, A1/B3 no
%   carrier, B2 carrier added after detection) have no analytic SNR model; they
%   are simulated at waveform level in oisl_kk_txcarrier.m, Section 8.
%
%   === INPUT PARAMETERS (all optional, with defaults) ===
%
%   Geometry / Channel:
%     w0        - Beam waist at Tx               [m]       (default: 0.1)
%     wd        - Detector aperture radius        [m]       (default: 0.07)
%     f         - Laser frequency                 [Hz]      (default: c/1530e-9)
%     hPL       - Terminal / path efficiency       [-]       (default: 0.9)
%     PT        - Total transmit power (carrier+signal) [W] (default: 0.3)
%     z         - Link distance                   [m]       (default: 1000e3)
%     sigma_s0  - Base pointing jitter            [m]       (default: 2.0)
%     k0        - Jitter growth rate               [-]       (default: 0.1)
%     d0        - Jitter reference distance        [m]       (default: 100e3)
%
%   Signal:
%     B         - Symbol rate                     [Hz]      (default: 10e9)
%     CSPR_dB   - Carrier-to-signal power ratio   [dB]      (default: 12)
%     beta      - RRC roll-off factor              [-]       (default: 0.05)
%     mod       - Modulation format                         (default: 'QPSK')
%                 Options: 'BPSK', 'OOK', 'QPSK'
%
%   Receiver:
%     eta       - Photodiode responsivity         [A/W]     (default: 0.5)
%     i_th      - Thermal noise current density   [A/sqrt(Hz)] (default: 20e-12)
%     RIN_dB    - Relative intensity noise        [dB/Hz]   (default: -150)
%     pre_on    - Pre-amplifier enabled            [bool]    (default: true)
%     G_dB      - Pre-amp gain                    [dB]      (default: 15)
%     NF_dB     - Pre-amp noise figure            [dB]      (default: 5)
%     Bf        - Optical filter bandwidth        [Hz]      (default: 30e9)
%     recon     - Add the heuristic KK reconstruction penalty [bool] (default: false).
%                 Not validated: 35-50 dB more pessimistic than the waveform
%                 simulation. The model without it is verified for CSPR >= 12 dB.
%
%   Simulation:
%     config    - Carrier/LO configuration; only 'tx_carrier' is modelled
%                 analytically (any other value raises an error)  (default: 'tx_carrier')
%     y_override - Pointing fraction y at which SNR / BER are reported [-]
%                  (default: [] = on-axis A0). Rates are always fading-averaged.
%
%   Unknown parameter names raise an error (a typo such as 'Pt' would otherwise
%   silently leave the default in place).
%
%   === OUTPUT STRUCT ===
%
%   result.params      - Copy of all parameters used
%   result.geom        - Beam geometry: lam, theta, wz, A0, sig, alpha
%   result.Prx_dBm     - Total received power at y=A0 [dBm]
%   result.Pc_dBm      - Carrier received power [dBm]
%   result.Ps_dBm      - Signal received power [dBm]
%   result.kappa       - CSPR (linear)
%   result.zeta        - KK noise enhancement factor (1+2/kappa)
%   result.PAPR_dB     - Signal PAPR for chosen modulation, 10^-4 CCDF point [dB]
%   result.snr_dB      - Total SNR at operating point [dB]
%   result.noise_budget - Struct: each field = individual noise SNR [dB]
%   result.rate_hard   - Ergodic rate, hard decision [bps]
%   result.rate_soft   - Ergodic rate, soft decision [bps]
%   result.rate_out10  - 10% outage rate [bps]
%   result.sens_dBm    - Receiver sensitivity [dBm] (Inf if KK fails)
%   result.sens_W      - Receiver sensitivity [W]
%   result.p_below_sens - Fraction of time the received power is below sensitivity
%   result.PT_min_1e2, result.PT_min_1e3 - Transmit power for 1e-2 / 1e-3 outage [W]
%   result.ber_theory  - Uncoded BER at the operating point y
%   result.gfec_dB     - SNR needed for the 3.8e-3 HD-FEC threshold [dB]
%   result.Prx_W, .snr, .y, .mod, .config - operating point echoes
%
%   Note: the noise model must stay identical to rx_snr() in oisl_kk_txcarrier.m.
%
%   === OPTISYSTEM INTEGRATION ===
%
%   This function is designed for OptiSystem MATLAB co-simulation.
%   OptiSystem provides realistic signal-level simulation; this function
%   provides the analytical channel model for comparison and validation.
%
%   Example usage from OptiSystem MATLAB Component:
%
%     % Get analytical predictions for current link configuration
%     r = oisl_kk_channel('PT', 0.3, 'z', 1500e3, 'mod', 'QPSK');
%     fprintf('Predicted SNR: %.1f dB, Rate: %.2f Gbps\n', ...
%             r.snr_dB, r.rate_hard/1e9);
%
%     % Compare OptiSystem BER with analytical
%     r = oisl_kk_channel('PT', PT_optisystem, 'f', f_optisystem);
%     analytical_ber = r.ber_theory;
%
%     % Sweep transmit power
%     for PT = 0.1:0.05:0.5
%         r = oisl_kk_channel('PT', PT, 'mod', 'QPSK');
%         fprintf('PT=%.2f W -> %.2f Gbps\n', PT, r.rate_hard/1e9);
%     end
%
%   Reference: ECE F266 Study Report, Chapter 3
%   Based on: Mecozzi et al., "Kramers-Kronig coherent receiver," Optica 2016
%
%   Date:   October 2026

    %% ---- Parse input ----
    P = default_params();
    if nargin == 1 && isstruct(varargin{1})
        fn = fieldnames(varargin{1});
        for i = 1:numel(fn), P = set_param(P, fn{i}, varargin{1}.(fn{i})); end
    elseif nargin > 0
        if rem(nargin, 2) ~= 0         % rem(), not mod(): mod is a variable below
            error('oisl_kk_channel:args', 'Parameters must be given as name-value pairs.');
        end
        for i = 1:2:nargin
            P = set_param(P, varargin{i}, varargin{i+1});
        end
    end
    if ~any(strcmp(P.mod, {'BPSK','OOK','QPSK'}))
        error('oisl_kk_channel:mod', 'mod must be ''BPSK'', ''OOK'' or ''QPSK'' (got ''%s'').', P.mod);
    end
    if ~strcmp(P.config, 'tx_carrier')
        error('oisl_kk_channel:config', ['Only config = ''tx_carrier'' has an analytic model. ' ...
            'Simulate ''%s'' with kk_link() in oisl_kk_txcarrier.m (Section 8).'], P.config);
    end

    % Extract simulation-only fields
    mod    = P.mod;
    config = P.config;
    y_ovr  = P.y_override;

    %% ---- Compute PAPR ----
    P.PAPR.BPSK = signal_papr_('BPSK', P.B, P.beta);
    P.PAPR.OOK  = signal_papr_('OOK',  P.B, P.beta);
    P.PAPR.QPSK = signal_papr_('QPSK', P.B, P.beta);

    %% ---- Beam geometry ----
    G = geom_(P, P.f, P.z);

    %% ---- Operating point ----
    k = 10^(P.CSPR_dB/10);
    if isempty(y_ovr)
        y = G.A0;           % on-axis (best case)
    else
        y = y_ovr;
    end
    Prx = P.hPL * P.PT * y;

    %% ---- SNR and noise budget ----
    [g, T] = rx_snr_(P, P.f, Prx, mod);

    % Build noise budget in dB
    fn = fieldnames(T); nb = struct();
    for i = 1:numel(fn)
        nb.(fn{i}) = -10*log10(T.(fn{i}));
    end

    %% ---- Rates ----
    rate_hard = R_erg_(P, P.f, P.PT, mod, 'hard');
    rate_soft = R_erg_(P, P.f, P.PT, mod, 'soft');
    rate_out  = R_out_(P, P.f, P.PT, mod, 0.1);

    %% ---- Sensitivity ----
    [sp, gfec] = sensitivity_(P, P.f, mod);

    %% ---- BER at operating point ----
    ber = ber_theory_(mod, g);

    %% ---- Assemble output ----
    result.params       = P;
    result.geom         = G;
    result.Prx_W        = Prx;
    result.Prx_dBm      = 10*log10(Prx/1e-3);
    result.Pc_dBm       = 10*log10(k/(1+k)*Prx/1e-3);
    result.Ps_dBm       = 10*log10(Prx/(1+k)/1e-3);
    result.kappa        = k;
    result.zeta         = 1 + 2/k;
    result.PAPR_dB      = 10*log10(P.PAPR.(mod));
    result.snr          = g;
    result.snr_dB       = 10*log10(g);
    result.noise_budget = nb;
    result.rate_hard    = rate_hard;
    result.rate_soft    = rate_soft;
    result.rate_out10   = rate_out;
    result.sens_dBm     = 10*log10(sp/1e-3);
    result.sens_W       = sp;
    result.p_below_sens = min(1, (sp/(P.hPL*P.PT*G.A0))^G.alpha);      % fraction of time below sensitivity
    result.PT_min_1e2   = sp/(P.hPL*G.A0*1e-2^(1/G.alpha));            % P_T for 1e-2 outage [W]
    result.PT_min_1e3   = sp/(P.hPL*G.A0*1e-3^(1/G.alpha));            % P_T for 1e-3 outage [W]
    result.ber_theory   = ber;
    result.gfec_dB      = gfec;
    result.config       = config;
    result.mod          = mod;
    result.y            = y;
end


function P = set_param(P, name, val)
% Only known parameter names are accepted (case-sensitive)
    if ~ischar(name) || ~isfield(P, name)
        error('oisl_kk_channel:param', 'Unknown parameter ''%s''. Valid names: %s', ...
            char(name), strjoin(fieldnames(P).', ', '));
    end
    P.(name) = val;
end


%% =========================================================================
%  DEFAULT PARAMETERS
%  =========================================================================
function P = default_params()
    c = 299792458;
    P.w0       = 0.1;           % beam waist [m]
    P.wd       = 0.07;          % detector aperture radius [m]
    P.f        = c / 1530e-9;   % laser frequency [Hz] (ITU-T C-band, 1530 nm)
    P.hPL      = 0.9;           % terminal / path efficiency
    P.PT       = 0.3;           % total transmit power [W]
    P.z        = 1000e3;        % link distance [m]
    P.sigma_s0 = 2.0;           % base pointing jitter [m]
    P.k0       = 0.1;           % jitter growth rate
    P.d0       = 100e3;         % jitter reference distance [m]
    P.B        = 10e9;          % symbol rate [Hz]
    P.CSPR_dB  = 12.0;          % CSPR [dB]
    P.beta     = 0.05;          % RRC roll-off
    P.eta      = 0.5;           % responsivity [A/W]
    P.i_th     = 20e-12;        % thermal noise [A/sqrt(Hz)]
    P.RIN_dB   = -150;          % RIN [dB/Hz]
    P.pre_on   = true;          % pre-amp on
    P.G_dB     = 15;            % pre-amp gain [dB]
    P.NF_dB    = 5;             % noise figure [dB]
    P.Bf       = 30e9;          % optical filter BW [Hz]
    P.pth      = 0;             % tracking threshold [W]
    P.recon    = false;         % heuristic reconstruction penalty (off: not validated)
    P.c        = c;
    P.q        = 1.602176634e-19;
    P.h        = 6.62607015e-34;
    % simulation fields
    P.mod       = 'QPSK';
    P.config    = 'tx_carrier';
    P.y_override = [];
end


%% =========================================================================
%  BEAM GEOMETRY
%  =========================================================================
function G = geom_(P, f, z)
    G.lam   = P.c / f;
    G.theta = G.lam / (pi * P.w0);
    G.wz    = z * tan(G.theta);
    G.A0    = 2 * P.wd^2 / G.wz^2;
    G.sig   = P.sigma_s0 * exp(P.k0 * z / P.d0);
    G.alpha = G.wz^2 / (4 * G.sig^2);
end


%% =========================================================================
%  RECEIVER SNR (heuristic reconstruction penalty only if P.recon)
%  =========================================================================
function [g, T] = rx_snr_(P, f, Prx, mod)
    k    = 10^(P.CSPR_dB/10);
    zeta = 1 + 2/k;
    hnu  = P.h * f;
    if P.pre_on
        Gr = 10^(P.G_dB/10); F = 10^(P.NF_dB/10);
        S  = hnu * (F*Gr - 1) / 2;
    else
        Gr = 1; S = 0;
    end
    Pc  = Gr * k/(1+k) * Prx;
    Ps  = Gr / (1+k) * Prx;
    sig = 2 * P.eta^2 * Pc .* Ps;
    Ibar = P.eta * (Pc + Ps + 2*S*P.Bf);
    RIN  = 10^(P.RIN_dB/10);

    T.shot    = zeta * (2*P.q*Ibar) .* P.B ./ sig;
    T.thermal = zeta * P.i_th^2 * P.B ./ sig;
    T.RIN     = zeta * RIN * (P.eta*(Pc+Ps)).^2 .* P.B ./ sig;
    if P.pre_on
        T.ASE     = (1+zeta) * S * P.B ./ Ps;
        T.ASExASE = zeta * 4 * P.eta^2 * S^2 * (P.Bf - P.B/2) * P.B ./ sig;
    end
    if nargin >= 4 && ~isempty(mod) && P.recon
        papr = P.PAPR.(mod);
        switch mod
            case 'QPSK',          pviol = exp(-k);
            case {'BPSK','OOK'},  pviol = erfc(sqrt(k/2));
        end
        T.recon = papr * max(pviol, 1e-30) + zeros(size(Prx));
    end
    fn = fieldnames(T); tot = 0;
    for i = 1:numel(fn), tot = tot + T.(fn{i}); end
    g = 1 ./ tot;
end


%% =========================================================================
%  MUTUAL INFORMATION
%  =========================================================================
function m = mi_(mod, g, kind)
    qf = @(x) 0.5*erfc(x/sqrt(2));
    h2 = @(p) -p.*log2(max(p,realmin)) - (1-p).*log2(max(1-p,realmin));
    if strcmp(kind, 'hard')
        switch mod
            case 'BPSK', m = 1 - h2(qf(sqrt(2*g)));
            case 'OOK',  m = 1 - h2(qf(sqrt(g)));
            case 'QPSK', m = 2*(1 - h2(qf(sqrt(g))));
        end
    else
        switch mod
            case 'BPSK', m = mi_bpsk_soft_(g);
            case 'OOK',  m = mi_bpsk_soft_(g/2);
            case 'QPSK', m = 2*mi_bpsk_soft_(g/2);
        end
    end
end


function I = mi_bpsk_soft_(g)
    persistent x w
    if isempty(x)
        n = 80; b = sqrt((1:n-1)/2); J = diag(b,1) + diag(b,-1);
        [V, D] = eig(J); x = diag(D); w = sqrt(pi)*V(1,:)'.^2;
    end
    I = zeros(size(g));
    for i = 1:numel(g)
        if g(i) <= 0, I(i) = 0; continue; end
        s2 = 1/(2*g(i)); yv = 1 + sqrt(2*s2)*x;
        zz = -2*yv/s2;
        l2 = (max(zz,0) + log1p(exp(-abs(zz)))) / log(2);
        I(i) = 1 - sum(w.*l2) / sqrt(pi);
    end
end


%% =========================================================================
%  THEORETICAL BER
%  =========================================================================
function p = ber_theory_(mod, g)
    qf = @(x) 0.5*erfc(x/sqrt(2));
    switch mod
        case 'BPSK', p = qf(sqrt(2*g));
        otherwise,   p = qf(sqrt(g));
    end
end


%% =========================================================================
%  ERGODIC RATE  (eq. 3.20)
%  =========================================================================
function R = R_erg_(P, f, PT, mod, kind)
    G = geom_(P, f, P.z);
    ylo = max(P.pth/(P.hPL*PT), G.A0*1e-30);
    if ylo >= G.A0, R = 0; return; end
    lA0 = log(G.A0);
    fun = @(t) P.B * mi_(mod, rx_snr_(P, f, P.hPL*PT*exp(t), mod), kind) ...
              .* (G.alpha * exp(G.alpha*(t - lA0)));
    R = integral(fun, log(ylo), lA0, 'RelTol', 1e-6, 'AbsTol', 1e-3);
end


%% =========================================================================
%  OUTAGE RATE  (eq. 3.21)
%  =========================================================================
function R = R_out_(P, f, PT, mod, eps_)
    G = geom_(P, f, P.z);
    y = G.A0 * eps_^(1/G.alpha);
    if y < P.pth/(P.hPL*PT), R = 0; return; end
    R = P.B * mi_(mod, rx_snr_(P, f, P.hPL*PT*y, mod), 'hard');
end


%% =========================================================================
%  SENSITIVITY  (BER = 3.8e-3 HD-FEC threshold)
%  =========================================================================
function [Ps, gfec] = sensitivity_(P, f, mod)
    gfec = fzero(@(gd) ber_theory_(mod, 10^(gd/10)) - 3.8e-3, [-10 30]);
    try
        Ps = 10^fzero(@(lp) 10*log10(rx_snr_(P, f, 10^lp, mod)) - gfec, [-14 0]);
    catch
        Ps = inf;
    end
end


%% =========================================================================
%  SIGNAL PAPR  (from RRC pulse-shaped waveform)
%  =========================================================================
function papr = signal_papr_(mod, B, beta)
% 10^-4 point of the instantaneous-power CCDF (2^17 symbols, fixed seed).
% Same definition as signal_papr() in oisl_kk_txcarrier.m.
    Nsym = 2^17; Nos = 8; N = Nsym*Nos; fs = Nos*B;
    fr = ((0:N-1).' - N*((0:N-1).' >= N/2)) * fs/N;
    switch mod
        case 'BPSK', pts = [1; -1];
        case 'OOK',  pts = [1; 0];
        case 'QPSK', pts = [1+1i; -1+1i; -1-1i; 1-1i]/sqrt(2);
    end
    rs = rng; rng(42);                 % fixed sequence, caller's RNG state restored below
    up = zeros(N,1); up(1:Nos:end) = pts(randi(numel(pts), Nsym, 1));
    rng(rs);
    s  = ifft(fft(up) .* rrc_(fr/B, beta));
    p  = sort(abs(s).^2 / mean(abs(s).^2));
    papr = p(ceil((1 - 1e-4)*N));
end


%% =========================================================================
%  ROOT RAISED COSINE
%  =========================================================================
function H = rrc_(fn, beta)
    af = abs(fn); H = zeros(size(fn));
    f1 = (1-beta)/2; f2 = (1+beta)/2;
    H(af <= f1) = 1;
    m = af > f1 & af <= f2;
    H(m) = sqrt(0.5*(1 + cos(pi/beta*(af(m) - f1))));
end
