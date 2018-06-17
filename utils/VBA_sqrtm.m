function S = VBA_sqrtm (C, inv)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% S = VBA_sqrtm (C, inv)
% quick derivation of the square root of a matrix
%
% IN:
%   - C: the entry matrix
%   - inv: binary flag for inverse square root (inv = true). 
%       Default = false.
%
% OUT:
%   - S: the (inverse) square root of C.
%
% /////////////////////////////////////////////////////////////////////////

if nargin < 2
    inv = false;
end

% degenerate case
if all (C(:) == 0)
    assert (~ inv, ...
        'VBA:expectedFinite', ...
        '*** VBA_sqrtm: cannot invert infinite precision.');
    S = 0 * C;
    
% diagonal case
elseif isequal (C, diag (diag (C)))
    dC = diag (C);
    if inv
        dC = 1 ./ dC;
    end
    S = diag (sqrt (dC));
    
% symmetric case
elseif issymmetric (C)
    assert (~ VBA_isWeird (C), ...
        'VBA:expectedFinite', ...
        '*** VBA_sqrtm: ill defined input.');
    [U, s] = svd (full (C));
    ds = diag (s);
    if inv
        ds = 1 ./ ds;
    end
    S = U * diag (sqrt (ds)) * U';
    
% general case
else
    error ('VBA_sqrtm: cannot compute sqrt of asymmetric matrix. See sqrtm().');
end
