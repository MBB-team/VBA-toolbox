function [a, b] = VBA_Create_NoisePrior(mu, sig)
%[ a,b ] = Create_NoisePrior( mu,sig )
%
% Converts mu and sig defining a prior distribution over noise standardard
% deviation s, to a gamma distriribution over the associated precision p
% defined by a and b with p = 1/s^2

% IN:
%   mu,sig: mean and standard deviation of prior distribution over s
% OUT:
%   a,b: Can either be a_sigma/b_sigma or a_alpha/b_alpha defining
%   prior distributions over measurement/system noise precision
% 
%  M. Eichenlaub 04/12/2019

if sig == 0     % Return Inf and 0 if stochstic inversion shall be disabled
    a = Inf;
    b = 0;
else

% Check mu and sig to prevent errors in numerical approximation        
assert(nargin==2,' *** VBA_Create_NoisePrior: Both mu and sigma must be provided')
assert(mu>5e-3,' *** VBA_Create_NoisePrior: mu must be greater than 5e-3')
assert(sig/mu>0.005,' *** VBA_Create_NoisePrior: sig must be greater than 0.005*mu')
assert(sig/mu<5,' *** VBA_Create_NoisePrior: sig must be smaller than 5*mu')

% Upper bound for a
a0 = 1/8*(1+sqrt(49+mu^4/sig^4+50*mu^2/sig^2)+mu^2/sig^2);

% Find a by minimizing the function D for 1<a<a0
a = fminbnd(@(x)D(x,mu,sig),1,a0);

% Calculate b
b = (mu./(exp(gammaln(a-0.5)-gammaln(a)))).^2;

% Obsolete because errors should be prevented earlier
% Check results, reject if error on either mu or sig is greater than 1%
% [m,s] = VBA_Convert_ab(a,b);
% if abs(m-mu)/mu>1e-2 || abs(s-sig)/mu>1e-2
%     error(' *** VBA_Create_NoisePrior: a and b could not be caluclated');
% end

end

end

function [ D ] = D(a,mu,sig)
% Function has to be minimized with respect to a

S = exp(2*(gammaln(a-0.5)-gammaln(a)));
D = mu^2./S - sig^2./((1./(a-1))-S);
D = log(D.^2+1);

end