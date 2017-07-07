function [sx,dsdx] = sparseTransform(x,P)
% operates a 'sparsify' mapping
% function [sx,dsdx] = sparseTransform(x,P)
% The sparse transform is applied to the input x, and is parameterized by
% P, where P is a smoothness parameter (default is 1).
sig = 1./(1+exp(-P*x));
sx = (x.^2).*(2.*sig-1);
dsdx = [];
if nargout==2
    dsigdx = P.*sig.*(1-sig);
    dsdx = 2*(x.^2).*dsigdx + 2*x.*(2.*sig-1);
    dsdx(x==0) = 1e-2; % for numerical purposes
    dsdx = diag(dsdx);
end

