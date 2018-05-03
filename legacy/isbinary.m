function flag = isbinary (X)
% legacy code
s = warning ('on');
warning ('*** The function `isbinary` is now deprecated and has been renamed `VBA_isBinary` (same syntax).') 
warning (s);

% fallback
flag = VBA_isBinary (X);