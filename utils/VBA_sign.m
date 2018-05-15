function [sx, dsdx] = VBA_sign (x, smoothness)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [sx, dsdx] = VBA_sign (x)
% smooth (ie. derivable) approximation of the sign() function
%
% IN:
%   - x: real value
%   - smoothness: optional parameter that specify the speed of transition
%     between -1 and 1. By default (1e-2), the change will occur for 
%     values of x between - 0.05 and 0.05. If Inf, then sx = sign (x) 
%
% OUT:
%   - sx: -1 if x is negative, +1 if it's positive. For values of x around
%     0, sx will progressively switch from -1 to 1.
%   - dsdx: derivative wrt x
%
% /////////////////////////////////////////////////////////////////////////

% parameters
% =========================================================================
% higher values means slower change from -1 to 1. If Inf, then sx = sign(x) 
if nargin < 2 || ~ isscalar (smoothness)
    smoothness = 1e-1;
end

% truncature to enforce non zero derivative
epsilon = 1e-9;

% sign approximation
% =========================================================================
% sigmoid scaled between -1 and 1
sx = 2 ./ (1 + exp(- x / smoothness)) - 1;

% derivative
% =========================================================================
if nargout < 2
    return
end

dsdx =  (sx + 1) .* (1 - 0.5 * (sx + 1)) / smoothness;

% stay positive!
dsdx(dsdx < epsilon) = epsilon;


