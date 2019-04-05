function [sx] = VBA_softmax (x, logTemperature, dim)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [sx] = VBA_softmax (x, logTemperature, dim)
% softmax function (generalized logistic mapping)
%
% IN:
%   - x: array of values to transform
%   - logTemperature: optional parameter (default = log(1)) that scales x
%       before the transformation (see below)
%   - dim: if x is multidimensional, dimension along which the softmax is 
%       applied (default to first non singleton dimension)
%
% OUT:
%   - sx: softmax transformation of x
%
% Background:
% ~~~~~~~~~~~ 
% 
%       softmax(x) = exp(x/tau) / sum_i [ exp(x_i/tau) ]
%
% with tau = exp(logTemperature) and i taken along dim
%
% /////////////////////////////////////////////////////////////////////////

% check inputs
% =========================================================================
if nargin < 2
    logTemperature = 0;
end

if nargin < 3
    % default to first non singleton dimension
    dim = find (size (x) > 1, 1, 'first');
else
    % check provided dimension is correct
    assert (size (x, dim) > 1, 'VBA_softmax: dim should be a nonsingleton dimension of x');
end

% compute softmax
% =========================================================================
% apply scaling
x = x / exp (logTemperature); 

% prevents numerical overflow
x = bsxfun (@minus, x, max (x, [], dim));

% actual mapping
ex = exp (x);
sx = bsxfun(@rdivide, ex, sum (ex, dim));








