function [lev] = lev_GLM(y,X)
% computes the log-evidence of a GLM (frequentist limit)
%  [lev] = lev_GLM(y,X)
% IN:
%   - y: dXn data array
%   - X: dXr design matrix
% OUT:
%   - lev: nX1 log-evidence array

try;X;catch;X=ones(size(y,1),1);end % default = mean
if isempty(X)
    d = size(y,1);
    r = 0;
    ldx = 0;
    yhat = 0;
else
    [d,r] = size(X);
    XtX = X'*X;
    ldx = VBA_logDet(XtX);
    yhat = X*VBA_inv(XtX)*X'*y;
end
e = y - yhat +eps;
n = size(y,2);
lev = zeros(n,1);
for i=1:n
    lev(i) = 0.5*(r-d)*log(2*pi) - 0.5*ldx + gammaln(0.5*(d-r)) + 0.5*(r-d)*log(0.5*e(:,i)'*e(:,i));
end

