function Y = VBA_quantile(X,p)
% Y = VBA_QUANTILE(X,P) returns quantiles of the values in the vector X.  
% VBA alternative to the quantile function from the statistics toolbox.
% P is a scalar or a vector of cumulative probability values. 
% Y is the same size as P, and Y(i) contains the P(i)-th quantile.


% check inputs
if isempty(X)
    Y = nan(numel(p));
end

X(isnan(X)) = [];
X = vec(X);

% compute quantiels
n = numel(X);
q = (0.5:1:(n-0.5))/n;
Y = interp1(q,sort(X),p);

Y(p==0) = min(X);
Y(p==1) = max(X);

