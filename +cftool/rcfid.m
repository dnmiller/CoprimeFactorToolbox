function [N, M] = rcfid(G)
% nrcf: Construct a right-coprime factorization with inner denominator
%
%   [N, M] = rcfid(G) 
% 
%       Construct a right-coprime factorization with inner denominator of the LTI system
%       G such that G = N/M and M'*M = I.

    assert(~isstable(G), 'System must be unstable');
    [A, B, C, D] = ssdata(balreal(minreal(G)));

    if isdt(G)
        [N, M] = drcfid(A, B, C, D, G.ts);
    elseif isct(G)
        [N, M] = crcfid(A, B, C, D);
    end
end

function [N, M] = crcfid(A, B, C, D)
% Compute a right coprime factorization with inner denominator for a
% continuous-time system.
% 
% Source: Zhou, et. al, Robust and Optimal Control, 1995.
    nu = size(B, 2);
    H = [A, -B * B'; zeros(size(A, 1)), -A'];
    X = gcare(H);
    F = -B' * X;
    
    fA = A + B * F;
    fB = B;
    fC = [F; C + D * F];
    fD = [eye(nu); D];
    
    CF = ss(fA, fB, fC, fD);
    
    M = CF(1:nu, :);
    N = CF(nu+1:end, :);
end

function [N, M] = drcfid(A, B, C, D, Ts)
% Compute a right coprime factorization with inner denominator for a
% discrete-time system.
% 
% Source: Chu, C.C., On Discrete Inner-Outer and Spectoral Factorizations,
% ACC, 1998.
    nu = size(B, 2);
    X = dare(A, B, zeros(size(A, 1)));
    R = eye(nu) + B' * X * B;
    F = -R \ B' * X * A;
    sqR = sqrtm(R);
    
    fA = A + B * F;
    fB = B / sqR;
    fC = [F; C + D * F];
    fD = [eye(nu); D] / sqR;
    
    CF = ss(fA, fB, fC, fD, Ts);
    
    M = CF(1:nu, :);
    N = CF(nu+1:end, :);
end