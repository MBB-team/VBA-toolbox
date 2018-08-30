function [x, dx] = VBA_abs (x, kappa)
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
% Note that at x = 0, VBA_abs (x) = 0.05 and dx = 1e-4;
%
% /////////////////////////////////////////////////////////////////////////

if nargin < 2
    kappa = 27.7259; % 0 => 0.05
end

% catch numericall unstable cases
if x > 5 
    x = abs(x);
    dx = 1;
    return
end

% compute derivative before we cange x...
if nargout > 1
    dx = tanh(0.5 * kappa * x) ;
    % ensure derivative is never flat
    if dx == 0
        dx = 1e-4;
    end
end

% approximate absolute value
x = - x + 2 * log1p (exp (kappa * x)) / kappa;
