function [C] = VBA_cov2corr (V)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [C] = cov2corr (V)
% derives correlation matrix from covariance matrix
%
% Let V(i,i), V(j,j) and V(i,j) be the variance of Xi, the variance of Xj
% and the covariance between Xi and Xj, respectively. Then the correlation
% between Xi and Xj is C(i,j) = V(i,j)./sqrt(V(i,i)*V(j,j)).
%
% IN:
%   - V: nxn covariance matrix
%
% OUT:
%   - C: nxn correlation matrix
%
% /////////////////////////////////////////////////////////////////////////

dv = diag (V);
C = V ./ sqrt (dv * dv');

