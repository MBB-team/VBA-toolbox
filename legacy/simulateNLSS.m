function [y,x,x0,eta,e,u] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0)
% legacy code
s = warning ('on');
warning ('*** The function `simulateNLSS` is now deprecated. Please use `VBA_simulate` instead (same syntax).') 
warning (s);

% fallback
try
    [y,x,x0,eta,e,u] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);
catch
    [y,x,x0,eta,e,u] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
end