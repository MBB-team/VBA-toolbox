function x = checkGX_binomial (x, lim)
% legacy code
s = warning ('on');
warning ('*** The function `checkGX_binomial` is now deprecated and has been renamed `VBA_finiteBinomial` (same syntax).') 
warning (s);

% fallback
x = VBA_finiteBinomial (x, lim);
