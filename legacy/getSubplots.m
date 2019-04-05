function getSubplots (ha, style)
% legacy code
s = warning ('on');
warning ('*** The function `getSubplots` is now deprecated and has been renamed `VBA_getSubplots` (same syntax).') 
warning (s);

% fallback
VBA_getSubplots (ha, style);