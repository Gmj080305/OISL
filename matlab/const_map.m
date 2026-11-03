function [pts, lab] = const_map(mod)
% Constellation points (unit average energy for BPSK/QPSK) and Gray bit labels
    switch mod
        case 'BPSK', pts = [1; -1];            lab = [0; 1];
        case 'OOK',  pts = [1; 0];             lab = [0; 1];
        case 'QPSK', pts = [1+1i; -1+1i; -1-1i; 1-1i]/sqrt(2); lab = [0 0; 1 0; 1 1; 0 1];
    end
end
