function [py, gridy] = VBA_empiricalDensity (y, pr)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [py, gridy] = VBA_empiricalDensity (y, pr)
% numerical derivation of the empirical distribution of y
%
% IN:
%   - y: n x p data matrix. The histograms are derived along the columns of
%       y, i.e. each column is considered as a vector of samples
%   - pr: the precision of the smoothing (gaussian) kernel (default = 1).
%
% OUT:
%   - py: (n-1) x p matrix of estimated emprical histograms
%   - gridy: (n-1) x p matrix of grid over which the emprical histograms are
%       estimated.
%
% Note: the following lines of code plots the empirical histogram:
% > [py, gridy] = VBA_empiricalDensity (y);
% > plot (gridy, py);
%
% /////////////////////////////////////////////////////////////////////////

% check input parameters
% =========================================================================
if nargin < 2
    pr = 1;
end

% initializations
% =========================================================================
% dimensions
[n, p] = size(y);

% smoothing kernel
kernel = exp (- pr * (- n / 2 : n / 2) .^ 2 / n);
kernel = kernel(:) / sum (kernel);

% memory allocation
py = zeros(n - 1, p);
gridy = zeros(n - 1, p);

% approximate density
% =========================================================================

% loop over columns
for i = 1 : p
    % get strictly increasing  values
    [sy, ecdf] = unique (sort (y(:, i)));
    
    % catch dirac case 
    if length(sy) == 1
        % center grid 
        midIdx = floor (n / 2);
        % grid around singularity
        bracket = sort ([0.99, 1.01] * sy);
        gridy(1 : midIdx - 1, i) = bracket(1);
        gridy(midIdx, i) = sy;
        gridy(midIdx + 1 : end, i) = bracket(2);
        % density to unit
        py(midIdx, i) = 1;
        continue;
    end
    
    % estimation grid
    gridyi = linspace (sy(1), sy(end), n);
        
    % approximate density
    ecdf = interp1 (sy, ecdf / n, gridyi);
        
    % smoothing
    py(:, i) = conv (diff (ecdf), kernel, 'same');
    py(:, i) = py(:, i) / sum (py(:, i));

    % reshape grid
    gridy(:, i) = (gridyi(1:end-1) + gridyi(2:end))/2;    
end
    
