function [pdf,cdf] = VBA_binomial(x,n,p)
% PDF and CDF of binomial distribution
% function [pdf,cdf] = VBA_binomial(x,n,p)
% Note: CDF is P(X<x), i.e. it excludes x, whereas PDF is P(X=x)!
% IN:
%   - x: number of successes
%   - n: total number of draws
%   - p: probability of a success
% OUT:
%   - pdf: P(X=x), where X is a binomial random variable
%   - cdf: P(X<x), where X is a binomial random variable

fullpdf = zeros(max(x),1);
for i=0:max(x)
    fullpdf(i+1) = nchoosek(n,i).*(p.^i).*((1-p).^(n-i));
end
pdf = fullpdf(x+1);
fullcdf = cumsum([0;fullpdf]);
cdf = fullcdf(x+1);