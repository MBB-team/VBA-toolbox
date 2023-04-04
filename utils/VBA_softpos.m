function [y, dydx] = VBA_softpos (x, varargin)
% SAFEPOSE is a derivable approximation of max(x,0).
% this function can be used to transform a parameter into positive numbers.
% IN: 
%   - x: value(s) to transform
%   - k: steepness of the transformation (default = 15.85, ie softpos(0) = 0.05);

% default values
k = 13.85; 
inverse = false;

% process arguments
for i = 1 : nargin-1
    if isnumeric(varargin{i})
        assert (k > 0, '*** VBA_softpos: k must be positive.');
        k = varargin{i};
    elseif strcmp (varargin{i}, 'inverse')
        inverse = true;
    else
        error ('*** VBA_softpos: invalid argument in position %d', i+1);
    end
end

if inverse
    y = log (exp (k * x) - 1) / k;
    return;
end

% main output
y = log (1 + exp (k * x)) / k;
% correct for saturation when x is too large
y(isinf (y)) = x(isinf (y));
    
if nargout > 1
    dydx = exp (k * x) ./ (exp (k * x) + 1);
end
end

