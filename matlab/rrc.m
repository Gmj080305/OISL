function H = rrc(fn, beta)
% Root-raised-cosine amplitude response, fn = f/B
    af = abs(fn); H = zeros(size(fn));
    f1 = (1-beta)/2; f2 = (1+beta)/2;
    H(af <= f1) = 1;
    m = af > f1 & af <= f2;
    H(m) = sqrt(0.5*(1 + cos(pi/beta*(af(m) - f1))));
end
