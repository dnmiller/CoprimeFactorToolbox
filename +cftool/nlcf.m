function [N, M] = nlcf(G)
% nlcf: Construct a normalized left-coprime factorization
%
%   [N, M] = nlcf(G) 
% 
%       Construct a normalized left coprime factorization of the LTI system
%       G such that G = M\N and M*M' + N*N' = I.

    [A, B, C, D] = ssdata(balreal(ss(minreal(G))));

    if isct(G)
        [N, M] = cnlcf(A, B, C, D);
    elseif isdt(G)
        [N, M] = dnlcf(A, B, C, D, G.Ts);
    end
end


function [N, M] = cnlcf(A, B, C, D)
% Compute a normalized left coprime factorization for the continuous-time
% case.
%
% Source: Zhou, et. al, Robust and Optimal Control, 1995.
    [ny, nu] = size(D);
    
    R = eye(nu) + D' * D;
    S = eye(ny) + D * D';
    V = A - B * D' / S * C;
    
    H = [V', -C' / S * C; -B / R * B', -V];
    Y = gcare(H);

    L = -(B * D' + Y * C') / S;
    
    sqS = sqrtm(S);
    
    A1 = A + L * C;
    B1 = [L, B + L * D];
    C1 = sqS \ C;
    D1 = sqS \ [eye(ny), D];
    
    NCF = ss(A1, B1, C1, D1);
    M = NCF(:, 1:ny);
    N = NCF(:, ny+1:end);
end


function [N, M] = dnlcf(A, B, C, D, Ts)
% Compute a normalized left coprime factorization for the discrete-time
% case.
%
% Source: Bongers, P.M., Heuberger, P.S., Discrete Normalized Coprime
% Factorization, 1990.
    ny = size(D, 1);
    
    S = eye(ny) + D * D';
    
    dA = A';
    dB = C';
    dQ = B * B';
    dR = S;
    dS = B * D';
    
    X = dare(dA, dB, dQ, dR, dS);
    
    Z = S + C * X * C';

    L = -(A * X * C' + B * D') / Z;
    sqZ = sqrtm(Z);
    
    fA = A + L * C;
    fB = [L, B + L * D];
    fC = sqZ \ C;
    fD = sqZ \ [eye(ny), D];
    
    NCF = ss(fA, fB, fC, fD, Ts);
    M = NCF(:, 1:ny);
    N = NCF(:, ny+1:end);
end