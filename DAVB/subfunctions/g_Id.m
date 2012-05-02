function [gx,dG_dX,dG_dPhi,d2G_dXdPhi] = g_Id(Xt,Phi,ut,inG)
% Identity observation mapping (partially observable)

n = size(Xt,1);
try
    ind = inG.ind;
catch
    ind = 1:n;
end
try
    scale = inG.scale;
catch
    scale = 1;
end
G = eye(n);
G = scale*G(ind,:);

try
    G = kron(G,ones(inG.k,1));
end


gx = G*Xt;
dG_dX = G';

if size(Phi,1) > 0
    dG_dPhi = zeros(size(Phi,1),size(G,1));
else
    dG_dPhi = [];
end
d2G_dXdPhi = [];