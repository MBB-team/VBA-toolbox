function [sx,dsdx,dsdP] = sparsify(x,P)
% legacy code
s = warning ('on');
warning ('*** The function `sparsify` is now deprecated. Please see `VBA_sparsifyPrior` for an alternative.') 
warning (s);

% fallback
[sx,dsdx,dsdP]  = VBA_sparsifyPrior(x, P);