function [N, M] = nrcf(G)
% nrcf: Construct a normalized right-coprime factorization
%
%   [N, M] = nrcf(G) 
% 
%       Construct a normalized right coprime factorization of the LTI system
%       G such that G = N/M and M'*M + N'*N = I.

    [A, B, C, D] = ssdata(balreal(ss(minreal(G))));
    
    if isct(G)
        [N, M] = cnrcf(A, B, C, D);
    elseif isdt(G)
        [N, M] = dnrcf(A, B, C, D, G.Ts);
    end
end


function [N, M] = cnrcf(A, B, C, D)
% Compute a normalized right coprime factorization for the continuous-time
% case.
%
% Source: Zhou, et. al, Robust and Optimal Control, 1995.
    [ny, nu] = size(D);
    
    R = eye(nu) + D' * D;
    S = eye(ny) + D * D';
    
    V = A - B / R * D' * C;
    
    H = [V, -B / R * B'; -C' / S * C, -V'];
    
    X = gcare(H);
    
    F = -R \ (B' * X + D' * C);
    sqR = sqrtm(R);
    
    fA = A + B * F;
    fB = B / sqR;
    fC = [F; C + D * F];
    fD = [eye(nu); D] / sqR;
    
    NCF = ss(fA, fB, fC, fD);
    
    M = NCF(1:nu, :);
    N = NCF(nu+1:end, :);
end

    
function [N, M] = dnrcf(A, B, C, D, Ts)
% Compute a normalized right coprime factorization for the discrete-time
% case.
%
% Source: Bongers, P.M., Heuberger, P.S., Discrete Normalized Coprime
% Factorization, 1990.
    nu = size(D, 2);
    
    R = eye(nu) + D' * D;
    
    dA = A;
    dB = B;
    dQ = C' * C;
    dR = R;
    dS = C' * D;
    
    X = dare(dA, dB, dQ, dR, dS);
    
    Z = R + B' * X * B;

    F = ((A' * X * B + C' * D) / Z)';
    sqZ = sqrtm(Z);
    
    fA = A - B * F;
    fB = B / sqZ;
    fC = [-F; C - D * F];
    fD = [eye(nu); D] / sqZ;
    
    NCF = ss(fA, fB, fC, fD, Ts);
    
    M = NCF(1:nu, :);
    N = NCF(nu+1:end, :);
end