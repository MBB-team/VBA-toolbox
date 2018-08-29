function lse = VBA_logSumExp (z, dim)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% lse = VBA_logSumExp (z)
% log of a sum of exponential, ie. log (sum (exp (z))
%
% IN:
%   - z: values to aggregate
%   - dim: if x is multidimensional, dimension along which the LSE is 
%       applied (default to first non singleton dimension)
%
% OUT:
%   - lse: log sum of exponentiated z
%
% Background:
% ~~~~~~~~~~~
% LSE function is a smooth approximation to the maximum function, ie.:
%   max(z) <= LSE(z) <= max(z) + log(n)
% but it also appears it a lot a numerical applications performed in the
% log domain space. However, LSE can quickly encounter numerical overflow
% problems and need to be computed with some extra care.
%
% /////////////////////////////////////////////////////////////////////////

% check inputs
% =========================================================================
if nargin < 2
    % default to first non singleton dimension
    dim = find (size (z) > 1, 1, 'first');
else
    % check provided dimension is correct
    assert (size (z, dim) > 1, 'VBA_logSumExp: dim should be a nonsingleton dimension of z');
end

% compute LSE
% =========================================================================
mz = max (z, [], dim);
z = bsxfun (@minus, z, mz);
lse = bsxfun(@plus, log (sum( exp (z), dim)), mz);
