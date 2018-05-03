function [x] = VBA_finiteBinomial(x, lim)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% x = VBA_finiteBinomial(x, lim)
% finesse binomial probabilities to avoid inifinite precision (p(x) = 0 or 1)
%
% IN:
%   - x: vector of binomial probabilities
%   - lim: probability threshold (default = 1e-8)
%
% OUT:
%   - x: vector of corrected binomial probabilities
%
% /////////////////////////////////////////////////////////////////////////

% set default
try
    lim;
catch
    lim = 1e-8;
end

% bound probabilities
x(x <= lim) = lim;
x(x >= 1 - lim) = 1 - lim;