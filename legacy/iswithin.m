function [flag] = iswithin (X, bounds)
% legacy code
s = warning ('on');
warning ('*** The function `iswithin` is now deprecated and has been renamed `VBA_isInRange` (same syntax).') 
warning (s);

% fallback
flag = VBA_isInRange (X, bounds);
