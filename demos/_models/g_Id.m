function [gx,dG_dX,dG_dPhi] = g_Id(Xt,Phi,ut,inG)
% Identity observation mapping (partially observable)

n = size (Xt, 1);

if ~ isfield (inG, 'ind')
    inG.ind = 1:n;
end
if ~ isfield (inG, 'scale')
    inG.scale = 1;
end

G = eye (n);
G = inG.scale * G(inG.ind, :);

if isfield(inG, 'k')
    G = kron (G, ones(inG.k, 1));
end


gx = G * Xt;
dG_dX = G';

if size (Phi, 1) > 0
    dG_dPhi = zeros (size (Phi, 1), size(G, 1));
else
    dG_dPhi = [];
end
