function [X] = VBA_sampleFromArbitraryP (p, gridX ,N)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [X] = VBA_sampleFromArbitraryP (p, gridX ,N)
% sample from an arbitrary 1D probability distribution
%
% IN:
%   - p: vector indicating the density evaluated along the grid (vector)
%   - gridX: the grid over which the density is evaluated (same size as p)
%   - N: the number of samples to be sampled (default = 1)
%
% OUT:
%   - X: Nx1 vector of samples
%
% /////////////////////////////////////////////////////////////////////////

% check inputs
% =========================================================================
% ensure proper dimensions
assert(isvector(p) && isvector(gridX), 'VBA_sampleFromP: density and grid must be vectors.');
assert(numel(p) == numel(gridX),  'VBA_sampleFromP: density and grid must have the same size.');

% check if density is really a density
assert(abs(sum(p) - 1) <= 1e-15, 'VBA_sampleFromP: density (p) must sum to 1!');

% fill in defaults
if nargin < 3
    N = 1;
end

% catch improper density
if any (isnan (p))
    X = nan(N, 1);
    return
end

% sample
% =========================================================================
% cumulative density
pcdf = cumsum (p(:));

% loop
X = zeros (N, 1);
for i = 1 : N
        below = find (rand <= pcdf, 1, 'first');
        X(i) = gridX(below);
end