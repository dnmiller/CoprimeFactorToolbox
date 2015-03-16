function result = isIoEqual(G0, G1, tol)
% isIoEqual: Test if two LTI systems are input-output equal.
% 
%   result = isIoEqual(G0, G1, tol)
% 
%       Return true if the two LTI systems 'G0' and 'G1' are input-output
%       equal and false otherwise. "equality" is determined by verifying
%       that the differences between the systems' impulse responses is
%       within the tolerance given by 'tol' (default is sqrt(eps)).
% 
    validateattributes(G0, {'ss', 'zpk', 'tf'}, {});
    validateattributes(G1, {'ss', 'zpk', 'tf'}, {});
    
    if nargin < 3
        tol = sqrt(eps);
    end

    assert(all(size(G0.d) == size(G1.d)), ...
        'G0 and G1 do not have same input/output dimension');
    
    assert((isct(G0) && isct(G1)) || (isdt(G0) && isdt(G1)), ...
        'G0 and G1 must be both discrete or both continuous');
    
    if isdt(G0)
        assert(G0.Ts == G1.Ts, 'G0 and G1 must have same sampling time');
    end

    [h0, t0] = impulse(G0);
    [h1, t1] = impulse(G1);
    
    if length(t0) > length(t1)
        [h1, ~] = impulse(G1, t0);
    else
        [h0, ~] = impulse(G0, t1);
    end
    
    maxErr = max(max(max(abs((h0 - h1)))));
    
    result = maxErr <= tol;
end