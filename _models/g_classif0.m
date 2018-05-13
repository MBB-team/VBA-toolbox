function [gx,dgdx,dgdp] = g_classif0(x,P,u,in)
% logistic-like observation function for VBA_classification
% function [gx,dgdx,dgdp] = g_classif0(x,P,u,in)
% IN:
%   - x: [useless]
%   - P: classification weights
%   - u: [useless]
%   - in: contains design matrix (in.X) and sparsify flag (in.sparse)
% OUT:
%   - gx: E[y|P]
%   - dgdx: [useless]
%   - dgdp: gradient of E[y|P] wrt P.

if in.sparse
    [sP, dsdP] = VBA_sparsifyPrior(P);    
else
    sP = P;
end
gx = sss(in.X'*sP);
dgdx = [];
dgdp = diag(gx.*(1-gx))*in.X';
dgdp = dgdp';
if in.sparse % for exploiting the analytical gradients from g_GLM
    dgdp = dsdP*dgdp;
end

function sx = sss(x)
sx = 1./(1+exp(-x));
sx(sx < 1e-8) = 1e-8;
sx(sx > 1-1e-8) = 1-1e-8;

