function [lev] = lev_GLM(y,X)
% computes the log-evidence of a GLM (frequentist limit)
%  [lev] = lev_GLM(y,X)
% IN:
%   - y: dXn data array
%   - X: dXr design matrix
% OUT:
%   - lev: nX1 log-evidence array
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

[d,r] = size(X);
XtX = X'*X;
yhat = X*pinv(XtX)*X'*y;
e = y - yhat;
n = size(y,2);
lev = zeros(n,1);
for i=1:n
    lev(i) = 0.5*(r-d)*log(2*pi) ...
        - 0.5*VBA_logDet(XtX) ...
        + gammaln(0.5*(d-r)) ...
        + 0.5*(r-d)*log(0.5*e(:,i)'*e(:,i));
end

