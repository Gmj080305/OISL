function D = fe_calibrate(Iref, SNL_dB, SME_dB, tap_samp)
% Fix the distortion coefficients so that, at the reference photocurrent Iref,
%   S_NL = 10 log10(S / N_NL) = SNL_dB   (paper eq. 13)
%   S_ME = 10 log10(S / N_ME) = SME_dB   (paper eq. 14)
% N is the power NOT explained by a best affine fit (gain + offset) of the
% distorted output on the input, i.e. the true distortion. The coefficients are
% then held fixed for every pointing state, as for a physical device.
    sig = std(Iref); u = (Iref - mean(Iref))/sig;
    D.sig_ref = sig; D.tap = tap_samp;
    D.a2 = solve_sdr(@(c) u + c*u.^2, u, SNL_dB);
    D.b  = solve_sdr(@(c) u + c*(circshift(u,tap_samp) + circshift(u,2*tap_samp)), u, SME_dB);
    D.SNL_dB = sdr_db(u, u + D.a2*u.^2);
    D.SME_dB = sdr_db(u, u + D.b*(circshift(u,tap_samp) + circshift(u,2*tap_samp)));
end

function c = solve_sdr(fun, u, target_dB)
% Bisection on log10(c) for a coefficient giving the requested signal-to-distortion ratio
    if isinf(target_dB), c = 0; return; end
    lo = -5; hi = 2;
    if sdr_db(u, fun(10^hi)) > target_dB
        warning('S = %.1f dB not reachable; using the strongest distortion (%.1f dB).', ...
            target_dB, sdr_db(u, fun(10^hi)));
        c = 10^hi; return;
    end
    for it = 1:60
        mid = (lo + hi)/2;
        if sdr_db(u, fun(10^mid)) > target_dB, lo = mid; else, hi = mid; end
    end
    c = 10^((lo + hi)/2);
end

function s = sdr_db(x, y)
    A = [ones(numel(x),1) x(:)]; c = A\y(:); res = y(:) - A*c;
    s = 10*log10(var(c(2)*x(:))/var(res));
end
