function [gx,dG_dX,dG_dPhi,d2G_dXdPhi] = g_sqrtSig(Xt,Phi,ut,inG)
% partially observable sigmoid mapping

n = size(Xt,1);
if ~exist('inG','var')
    inG = [];
end
try
    ind1 = inG.ind;
catch
    ind1 = 1:n;
end
G = eye(n);
G = G(ind1,:);

[Sx,dsdx] = sigm(Xt,inG);

gx = G*Sx(:).^2;
dG_dX = 2*[G*diag(Sx(:).*dsdx(:))]';
dG_dPhi = [];
d2G_dXdPhi = [];

