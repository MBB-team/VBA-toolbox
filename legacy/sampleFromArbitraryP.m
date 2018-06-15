function [X] = sampleFromArbitraryP (p, gridX ,N)
% legacy code
s = warning ('on');
warning ('*** The function `sampleFromArbitraryP` is now deprecated. Please see `VBA_random` for an alternative.') 
warning (s);

% fallback
if nargin < 3
    N = 1;
end

if isvector (gridX)
    N = {N, 1};
else
    N = {N};
end
X = VBA_random ('Arbitrary', p, gridX, N{:});