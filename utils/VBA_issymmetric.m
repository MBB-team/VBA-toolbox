function [flag] = VBA_issymmetric (X)

% // VBA toolbox //////////////////////////////////////////////////////////
%
% [flag] = VBA_issymmetric (X)
% check if X is a symmetric matrix
%
% IN:
%   - X: matrix
% OUT:
%   - flag: - true if X is symmetric
%
% /////////////////////////////////////////////////////////////////////////

try
    flag = issymmetric (X);
catch
    flag = isequal (X, X');
end