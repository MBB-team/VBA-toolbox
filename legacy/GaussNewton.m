function [opt,sigma,out] = GaussNewton(fname,init,options)
% legacy code
s = warning ('on');
warning ('*** The function `GaussNewton` is now deprecated and has beend renamed `VBA_GaussNewton`.') 
warning (s);

% fallback
try
    options;
catch
    options = struct;
end
[opt,sigma,out] = VBA_GaussNewton(fname,init,options);