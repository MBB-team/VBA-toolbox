function [sx,dsdx,dsdP] = sparsify(x,P)
% operates a 'sparsify' mapping
% function [sx,dsdx] = sparsify(x,P)
% The sparse transform is applied to the input x, and is parameterized by
% P, where P is the sparsity parameter.
sig = 1./(1+exp(-x));
sx = (abs(x).^exp(P)).*(2.*sig-1);
if nargout>=2
    dsigdx = sig.*(1-sig);
    dsdx = 2*(abs(x).^exp(P)).*dsigdx + exp(P).*abs(x).^(exp(P)-1).*(2.*sig-1).^2;
    dsdx(x==0) = 1e-2; % for numerical purposes
    dsdx = diag(dsdx);
    dsdP = sx.*log(abs(x+1e-2)).*exp(P);
end

function gax = myAbs(x)
a = 1;
gax = -x + 2*log(1+exp(a.*x))./a;