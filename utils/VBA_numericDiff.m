function [dfdx, fx] = VBA_numericDiff(fName, idxArg2Diff, varargin)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [dfdx, fx] = numericDiff(fname, idxArg2Diff, fArg1, fArg2, ...)
% numerical evaluation of derivatives
%
% This function evaluates numerically the derivatives of the function
% 'fname' with respect to its 'idxArg2Diff' input argument at the point 
% defined by the input arguments [fArg1, fArg2, ..., fArgN]
%
% IN:
%   - fname: handle or name of the function to be differentiated.
%     ie either its name or a function handle. The latter can be very useful
%     if the function is itself a subfunction of a function!
%   - idxArg2Diff: the index of the argument to differentiate the function
%     with (1 <= idxArg2Diff <= numel(varargin))
%   - fArg1, fArg2, ..., fArgN: list of arguments which is required to call
%     the function 'fname', and which defines the ordinate at which the 
%     derivative will be numerically evaluated
%
% OUT:
%   - dfdx: an m x p matrix containing the numerical differentiation
%     of the function evaluated at {fArg1, ..., fArgN}, where p is the 
%     number of elements of the function's first output and m is the number
%     of elements of the idxArg2Diff input argument. If the function output
%     or the differentiated argument are in matrix form, they will be
%     vectorized first.
%   - f: the function evaluation at x
%
% NB: Mixed partials can be obtained by recursive call of this routine,
% e.g. consider:
%
% dfdx_idx_j = numericDiff(@numericDiff,i+2,fname,j,arg1,arg2,...,argn)
% dfdx_idx_j = reshape(dfdx_idx_j,mi,mj,p)
%
% The first line evaluates the numerical derivative of the function
% numericDiff(fname,j,arg1,arg2,...,argn) wrt to its (i+2)th entry, which
% is arg_i.
% The second line reshapes the output such that the first dimension is the
% dimension of argi (mi), the second is of the dimension of argj (mj), and
% the last dimension is the one of the function output itself (p).
%
% /////////////////////////////////////////////////////////////////////////

% perturbation scale
epsilon = 1e-4;

% shortcut
fArgs = varargin;

% evaluate function at the specified argument
% =========================================================================
try
    fx = VBA_vec(fName(fArgs{:}));
catch 
    message = sprintf(...
        '*** VBA_numericDiff: Can not call %s with the provided arguments (nArgs = %d).', ...
        func2str(fName), numel(fArgs));
    error(message);
end

% evaluate the perturbation df of the function fname in the neighbourhood 
% of the specified argument (ie at x + dx) and calculate the derivative wrt x
% =========================================================================

% pre-allocate the variables
m = numel (fArgs{idxArg2Diff});
p = numel (fx);
dfdx = zeros (m, p);

% loop over the dimension of idxArg2Diff-th argument
dx = epsilon * fArgs{idxArg2Diff};
dx (abs(dx) <= eps) = epsilon;

% compute effect of pertubations
for i = 1 : m
    xpdx = fArgs;
    xpdx{idxArg2Diff}(i) = xpdx{idxArg2Diff}(i) + dx(i);
    dfdx(i,:) = (VBA_vec (fName (xpdx{:})) - fx)' / dx(i);
end

