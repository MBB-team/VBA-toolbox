function [flag] = isweird (X)
% legacy code
s = warning ('on');
warning ('*** The function `isweird` is now deprecated and has been renamed `VBA_isWeird` (same syntax).') 
warning (s);

% fallback
flag = VBA_isWeird (X);
