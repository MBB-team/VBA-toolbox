function p = VBA_ncfcdf(x,nu1,nu2,delta)
% NCFCDF Noncentral F cumulative distribution function (cdf).
% function p = VBA_ncfcdf(x,nu1,nu2,delta)
% Returns the noncentral F cdf with numerator degrees of freedom (df), NU1,
% denominator df, NU2, and noncentrality parameter, DELTA, at the values in
% X. The size of P is the common size of the input arguments. A scalar
% input functions as a constant matrix of the same size as the other
% inputs.

n = numel(x);
J = 0:1e3;
summand = exp(-0.5*delta).*((0.5*delta).^J)./factorial(J);
summand(isnan(summand)) = 0; % cf. limit e^n/n!
p = NaN(size(x));
for i=1:n
    Ix = betainc(nu1*x(i)./(nu2+nu1*x(i)),0.5*nu1+J,0.5*nu2);
    Ix(isnan(Ix)) = 0;
    p(i) = sum(summand.*Ix);
end

