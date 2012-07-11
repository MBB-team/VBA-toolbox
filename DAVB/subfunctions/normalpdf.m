function [p] = normalpdf(x,E,V)
% computes the normal probability density function
% [p,E,V] = normalpdf(x,E,V)
% IN:
%   - x: the value of the random variable that is normal-distributed
%   - E/V: first- and second-order moments of the normal density
% OUT:
%   - p: the beta probability density function evaluated at x

[k,n] = size(x);

logp = zeros(1,n);
logp(1:n) = -(k/2).*log(2*pi) - 0.5*VBA_logDet(V);
iV = VB_inv(V);
for i=1:n
    e = x(:,i)-E;
    logp(i) = logp(i) -0.5*e'*iV*e;
end
p = exp(logp);