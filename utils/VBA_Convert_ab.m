function [mu, sig] = VBA_Convert_ab(a, b)
% [mu,sig] = Convert_ab(a,b)
%
% Converts a and b defining a gamma distributed precision p to mean (mu) and
% standard deviation (sig) of the distribution over the associated standard
% deviation s, i.e. s = 1/sqrt(p)
% 
% IN:
%   a,b: Can either be a_sigma/b_sigma or a_alpha/b_alpha defining
%   posterior distributions over measurement/system noise precision
% OUT:
%   mu,sig: mean and standard deviation of distribution over s
% 
% Be aware that this is only possible for a>1
%
%  M. Eichenlaub 04/12/2019

if isinf(a) % Return 0 if system noise was switched off, i.e. a_sigma = Inf
    mu = 0;
    sig = 0;
elseif a<=1
    error('*** VBA_Convert_ab: a must be greater or equal to 1');
else    
    mu = sqrt(b)*exp(gammaln(a-0.5)-gammaln(a));
    sig = sqrt(b/(a-1)-b*exp(2*(gammaln(a-0.5)-gammaln(a))));
end

 


