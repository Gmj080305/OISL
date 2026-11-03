function o = rx_metrics(L, Npilot)
% Pilot-aided one-tap gain/phase estimate (the first Npilot symbols are known
% pilots, as a real receiver would have). SNR, BER, EVM and phase error are then
% measured on the payload only, so no genie knowledge of the data is used.
    a = L.a; r = L.r; pts = L.pts; lab = L.lab;
    np = 0; if isfield(L, 'npost'), np = L.npost; end      % end pilots (frequency recovery)
    p  = 1:Npilot; d = Npilot+1:numel(a)-np;
    Gc = (a(p)'*r(p))/(a(p)'*a(p));
    req = r/Gc;                                            % equalised, unit-energy scale
    e   = req(d) - a(d);
    o.snr  = mean(abs(a(d)).^2)/mean(abs(e).^2);
    [~, dec] = min(abs(req(d) - pts.').^2, [], 2);
    o.nerr  = sum(sum(lab(L.idx(d),:) ~= lab(dec,:)));
    o.nbits = numel(lab(L.idx(d),:));
    o.ber   = o.nerr/o.nbits;
    o.evm   = sqrt(mean(abs(e).^2)/mean(abs(a(d)).^2));
    o.phase_err = angle(req(d)./a(d));
    o.sigma_phi = std(o.phase_err);
    o.sym = req(d); o.tx = a(d); o.idx = L.idx(d);
    o.sym_plot = o.sym/sqrt(mean(abs(o.sym).^2));          % unit-power view (same axes for all)
end
