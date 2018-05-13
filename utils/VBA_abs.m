function [x, dx] = VBA_abs (x)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [x] = VBA_abs (x)
% smooth (ie. derivable) approximation of the abs() function
%
% IN:
%   - x: real value
%
% OUT:
%   - x: absolute value of x (approximated)
%   - dx: derivative wrt x
%
% Note that at x = 0, VBA_abs (x) = 0.01 and dx = 1e-4;
%
% /////////////////////////////////////////////////////////////////////////

z = 138.6294;

% catch numericall unstable cases
if x > 5 
    x = abs(x);
    dx = 1;
    return
end

% compute derivative before we cange x...
if nargout > 1
    dx = tanh(0.5 * z * x) ;
    % ensure derivative is never flat
    if dx == 0
        dx = 1e-4;
    end
end

% approximate absolute value
x = - x + 2 * log1p (exp (z * x)) / z;
