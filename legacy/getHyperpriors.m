function [a, b] = getHyperpriors(y, p_min, p_max)
% legacy code
s = warning ('on');
warning ('*** The function `getHyperpriors` is now deprecated. Please see `VBA_guessHyperpriors` for an alternative.') 
warning (s);

% fallback
[a, b] = VBA_guessHyperpriors (y, [p_min, p_max]);