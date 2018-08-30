function [sx, dsdx, dsdp] = VBA_sparsifyPrior (x, logExponent, smoothness)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [sx, dsdx, dsdp] = VBA_sparsifyPrior (x, varargin)
% parameter transformation that emulates Laplace priors (L1-norm)
%
% IN:
%   - x: input value to be transformed
%   - logExponent: optional sparsity parameter. If P = 0, then the mapping 
%     is regular and no sparsity is emulated. If P=log(2), then the mapping
%     emulates L1-norm like priors.
%   - smoothness: optional smoothness parameter of the sign approximation 
%     (default = 1)
%
% In the general case, you won't have to change any of the optional
% parameters. Those options allows to explore the properties of the
% transformation and should not be changed.
%
% OUT:
%   - sx: transformed value
%   - dsdx: gradient of sx wrt x
%   - dsdP: gradient of sx wrt logExponent parameter
%
% /////////////////////////////////////////////////////////////////////////

% check parameters
% =========================================================================
if nargin < 2
    logExponent = log(2);
end

if nargin < 3
    smoothness = []; % use VBA_sign default
end

% shortcuts
% =========================================================================
[signX, d_signX] = VBA_sign(x, smoothness);
[absX, d_absX] = VBA_abs (x);
exponent = exp(logExponent);

% compute transformation
% =========================================================================

sx = signX .* (absX .^ exponent);

% derivatives
% =========================================================================

% wrt to x
% -------------------------------------------------------------------------
if nargout < 2
    return;
end

dsdx = d_signX .* (absX .^ exponent) ...
     + exponent .* absX .^ (exponent - 1) .* signX .* d_absX;      

% wrt to logExponent
% -------------------------------------------------------------------------
if nargout < 3
    return;
end

dsdp = (absX .^ exponent) .* signX .* log(absX) .* exponent;
