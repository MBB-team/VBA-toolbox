function [dfdx,f] = numericDiff(fname,nArg2diff,varargin)
% numerical evaluation of derivatives
% [dfdx] = numericDiff(fname,nArg2diff,arg1,arg2,...,argn)
% This function evaluates numerically the derivatives of the function
% called fname. 
% IN:
%   - fname: a variable which refers to the function to be differentiated,
%   ie either its name or a function handle. !The latter can be very useful
%   if the function is itself a subfunction of a function!
%   - nArg2diff: the index of the argument to differentiate the function
%   with (1 <= nArg2diff <=length(varargin))
%   - arg1,arg2,...,argn: the list of arguments which is required to call
%   the function fname (at which the derivative will be numerically
%   evaluated)
% OUT:
%   - dfdx: an mxp matrix containing the numerical differentiation
%   of the function at x=Arg, where p is the vectorized dimension of the
%   function fname, and m is the vectorized dimension of argument with
%   which the function is differentiated
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


Arg = varargin;
clear varargin

nArgs = length(Arg);
nArgIn = nargin(fname);

if nArgIn<0 || isequal(nArgs,nArgIn)

    %-- evaluate function fname at the specified argument (at x)  --%
    f = feval(fname,Arg{:});
    fx = vec(f);   
    
    %-- evaluate the perturbation df of the function fname in the
    % neighbourhood of the specified argument (ie at x + dx) and
    % calculate the derivative wrt x                              --%

    % pre-allocate the variables
    m = numel(Arg{nArg2diff});
    dfdx = zeros(m,size(fx,1));
    % loop over the dimension of x
    for i = 1:m
        xpdx = Arg;
        dx = 1e-4*xpdx{nArg2diff}(i);
        if abs(dx) <= eps    % whoops !!
            dx = 1e-4;
        end
        xpdx{nArg2diff}(i) = xpdx{nArg2diff}(i) + dx;
        dfdx(i,:) = (vec(feval(fname,xpdx{:})) - fx)'./dx;
    end
    
else
    
    error('numericDiff: wrong number of input arguments...')

end


function vx = vec(X)
vx = X(:);  % returns column vector

