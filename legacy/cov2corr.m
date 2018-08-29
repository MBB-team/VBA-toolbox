function y = cov2corr(x)
% legacy code
s = warning ('on');
warning ('*** The function `cov2corr` is now deprecated. Please see `VBA_cov2corr` (same syntax).') 
warning (s);

% fallback
y = VBA_cov2corr (x);