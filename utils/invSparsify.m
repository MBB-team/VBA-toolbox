function [isx] = invSparsify(x,P)
% operates the inverse 'sparsify' mapping (see sparsify.m)
isx = (abs(x).^exp(-P)).*sign(x);