function x = checkGX_binomial(x,lim)

% finesses 0/1 (inifinite precision) binomial probabilities

% function x = checkGX_binomial(x)
% IN:
%   - x: vector of binomial probabilities
%   - lim: probability threshold
% OUT:
%   - x: vector of corrected binomial probabilities
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

if nargin==1
    lim = 1e-8;
end
x(x<=lim) = lim;
x(x>=1-lim) = 1-lim;