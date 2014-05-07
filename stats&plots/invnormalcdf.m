function [x] = invnormalcdf(p,m,v)
% computes the inverse normal cumulative density function
% [p,E,V] = invnormalcdf(p,E,V)
% IN:
%   - p: the cumulative density value
%   - E/V: first- and second-order moments of the normal density
% OUT:
%   - x: the value of the random variable X such that P(X>x|E,V) = p

p = p(:);
n = length(p);

% calculate normal cdf on a grid
s = sqrt(v);
dx = s*1e-4;
gridx = -8*s:dx:8*s;
gridx = gridx + m;
pdf = exp(-0.5*(gridx-m).^2./v);
pdf = pdf./sum(pdf);
cdf = cumsum(pdf);

% find x such that P(X>x|E,V) = p
x = zeros(size(p));
for i=1:n
    [mp,indp] = min(abs(cdf-p(i)));
    x(i) = gridx(indp);
end