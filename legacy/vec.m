function vx = vec(X)
% legacy code
s = warning ('on');
warning ('*** The function `vec` is now deprecated and has beend renamed `VBA_vec`.') 
warning (s);

% fallback
vx = VBA_vec (X);