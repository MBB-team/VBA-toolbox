function [X] = sampleFromArbitraryP (p, gridX ,N)
% legacy code
s = warning ('on');
warning ('*** The function `sampleFromArbitraryP` is now deprecated and has beend renamed `VBA_sampleFromArbitraryP`.') 
warning (s);

% fallback
if nargin < 3
    N = 1;
end
[X] = VBA_sampleFromArbitraryP (p, gridX ,N);