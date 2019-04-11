function V = VBA_getVar(Sigma,T)
% gets the time-dependent variances
% function V = VBA_getVar(Sigma,T)
% This function gets the (possibly time-dependent) variance stored in a
% cell array of covariance matrices.
% IN:
%   - Sigma: the cell array of covariance matrices
%   - T: the index of the end of the time series
% OUT:
%   - V: a nXT matrix containing the variances in the leading diagonal of
%   the nXn matrices contained in the cell array Sigma.


if isempty (Sigma)
    V = [];
elseif isnumeric (Sigma)
    V = diag(Sigma);
else
    if nargin < 2
        T = length(Sigma);
    end
    n = size(Sigma{1}, 1);
    V = zeros(n,T);
    for t=1:T
        V(:,t) = diag(Sigma{t});
    end
end