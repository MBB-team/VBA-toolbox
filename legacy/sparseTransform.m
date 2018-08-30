function [sx,dsdx] = sparseTransform(x,P)
% legacy code
s = warning ('on');
warning ('*** The function `sparseTransform` is now deprecated. Please see `VBA_sparsifyPrior` for an alternative.') 
warning (s);

% fallback
[sx,dsdx] = VBA_sparsifyPrior(x, log(2), 1/P);