function Els = VBA_Elogsig (m, V)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% Els = VBA_Elogsig (m, V)
% semi-analytical approx to E[log(sig(x))] when x~N(m,V)

% This function computes a semi-analytical approximation to the expected
% log-sigmoid mapping of x, when x is a gaussian variable with mean m and
% variance V.

% IN:
%   - m: mean of x
%   - V: variance of x
%
% OUT:
%   -Els: E[log(sig(x))]
%
% See "Semi-analytical approximations to statistical moments of sigmoid and
% softmax mappings of normal variables" (Daunizeau 2017).
%
% /////////////////////////////////////////////////////////////////////////

% Constants
a = 0.205;
b = - 0.319;
c = 0.781;
d = 0.870;

% Computation
tmp = (m + b .* V .^ c) ./ sqrt (1 + a .* V .^ d);
Els = log(1 ./ (1 + exp (- tmp)));