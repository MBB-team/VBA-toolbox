function Els = Elogsig(m,V)
% computes semi-analytical approx to E[log(sig(x))] when x~N(m,V)
% function Els = Elogsig(m,V)
% This function computes a semi-analytical approximation to the expected
% lo-sigmoid mapping of x, when x is a gaussian variable with mean m and
% variance V.
% IN:
%   - m: mean of x
%   - V: variance of x
% OUT:
%   -Els: E[log(sig(x))]
% See "Semi-analytical approximations to statistical moments of sigmoid and
% softmax mappings of normal variables" (Daunizeau 2017).

a = 0.205;
b = -0.319;
c = 0.781;
d = 0.870;

tmp = (m+b.*V.^c)./sqrt(1+a.*V.^d);
Els = log(VBA_sigmoid(tmp));