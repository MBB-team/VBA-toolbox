function [H, pv] = VBA_kstest(x,alpha)
% derives Kolmogorov-Smirnov test of gaussianity
% [H, pv] = VBA_kstest(x)
% This was adapted from the function normalitytest.m, which is available
% from MATLAB's flie exchange.

try;alpha;catch,alpha=0.05;end % test size
x = x(:);
n = length(x);
y = sort(x);
y = bsxfun(@minus,y, mean(y));
sigma0 = std(y);
sigma0(sigma0==0) = 1;
y = bsxfun(@rdivide,y, sigma0);
empcdf = [1:n]'/n;
empcdf2 = ([1:n]-1)'/n;
fx = VBA_spm_Ncdf(y);
KSz = max([abs(fx-empcdf);abs(fx-empcdf2)]);
s = -20:1:20;
a = (-1).^s.*exp(-2*(s.*sqrt(n)*KSz).^2); 
pv = 1-sum(a);
H = (pv<alpha);
