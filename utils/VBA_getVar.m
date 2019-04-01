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


if iscell(Sigma)
    if nargin < 2
        T = length (Sigma);
    end
    V = cell2mat (cellfun (@diag, Sigma(1 : T), 'UniformOutput', false));
else
    V = diag(Sigma);
end