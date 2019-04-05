function [gx] = g_DG2(Xt,Phi,ut,inG)
% double-Gamma function HRF observation function

grid = inG.grid;

a1 = exp(Phi(1));
b1 = exp(Phi(2));
a2 = exp(Phi(3));
b2 = exp(Phi(4));
lambda1 = Phi(5);
lambda2 = Phi(6);
CST = Phi(7);

gx = zeros(size(grid,1),1);
gx = lambda1.*(b1.^a1./gamma(a1)).*grid.^(a1-1).*exp(-b1.*grid) ...
    + lambda2.*(b2.^a2./gamma(a2)).*grid.^(a2-1).*exp(-b2.*grid) ...
    + CST;