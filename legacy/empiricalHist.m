function [py, gridy] = empiricalHist (y, pr)
% legacy code
s = warning ('on');
warning ('*** The function `empiricalHist` is now deprecated and has beend renamed `VBA_empiricalDensity`.') 
warning (s);

% fallback
[py, gridy] = VBA_empiricalDensity (y, pr);