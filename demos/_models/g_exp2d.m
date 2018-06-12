function [gx] = g_exp2d(x,phi,u,in)

gkaf = in.gkaf;
gkas = in.gkas;
% 
% X = [gkaf,gkas,gkaf.^2,gkas.^2,gkaf.*gkas];
% X = [X,zeros(size(X,1),1)];
% 
% gx = X*phi;
% dgdx = [];
% dgdphi = X';


gx = phi(1).*exp(0.5.*gkaf).*exp(-abs(gkas-0.5.*gkaf));
% gx = real(gx);