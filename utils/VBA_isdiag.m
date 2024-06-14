function flag = VBA_isdiag (A)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% flag = VBA_isdiag (A)
% check if a matrix is diagonal
%
% IN:
%   - A: matrix
% OUT:
%   - flag: boolean flag, true if the matrix A is diagonal
%
% Notes:
% ~~~~~~
% this is a fallback for isdiag on Matlab < 2014a
%
% /////////////////////////////////////////////////////////////////////////

try
    flag = isdiag (A);
    return
catch
    flag = all (VBA_vec (A == diag (diag (A))));
end

