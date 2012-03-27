function [gx,dG_dX,dG_dPhi,d2G_dXdPhi] = g_DoubleGamma(Xt,Phi,ut,inG)
% double-Gamma function HRF observation function
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
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
    + lambda2.*(b2.^a2./gamma(a2)).*grid.^(a2-1).*exp(-b2.*grid)...
    +CST;

dG_dPhi = numericDiff(@g_DG2,2,Xt,Phi,ut,inG);


% build the observation function gradients
% 
% dG_dPhi(:,1) = lambda1.*b1.^a1.*log(b1)./gamma(a1).*grid.^(a1-1).*exp(-b1.*grid)...
%     -b1.^a1./gamma(a1).*grid.^(a1-1).*exp(-b1.*grid).*Psi(a1)...
%     +b1.^a1./gamma(a1).*grid.^(a1-1).*log(grid).*exp(-b1.*grid);
% dG_dPhi(:,2) = lambda1.*b1.^a1.*a1./b1./gamma(a1).*grid.^(a1-1).*exp(-b1.*grid)...
%     -b1.^a1./gamma(a1).*grid.^(a1-1).*grid.*exp(-b1.*grid);
% 
% dG_dPhi(:,3) = lambda2.*b2.^a2.*log(b2)./gamma(a2).*grid.^(a2-1).*exp(-b2.*grid)...
%     -b2.^a2./gamma(a2).*grid.^(a2-1).*exp(-b2.*grid).*Psi(a2)...
%     +b2.^a2./gamma(a2).*grid.^(a2-1).*log(grid).*exp(-b2.*grid);
% dG_dPhi(:,4) = lambda2.*b2.^a2.*a2./b2./gamma(a2).*grid.^(a2-1).*exp(-b2.*grid)...
%     -b2.^a2./gamma(a2).*grid.^(a2-1).*grid.*exp(-b2.*grid);
% 
% dG_dPhi(:,5) = (b1.^a1./gamma(a1)).*grid.^(a1-1).*exp(-b1.*grid);
% dG_dPhi(:,6) = (b2.^a2./gamma(a2)).*grid.^(a2-1).*exp(-b2.*grid);
% dG_dPhi = dG_dPhi';


dG_dX = [];
d2G_dXdPhi = [];






