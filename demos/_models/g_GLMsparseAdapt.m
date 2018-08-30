function [gx,dgdx,dgdP] = g_GLMsparseAdapt(x,P,u,in)


[sP, dsdx, dsdp] = VBA_sparsifyPrior (P(1:end-1), P(end));

[gx, ~, dgdP] = g_GLM (x, sP, u, in);

dgdx = [];

dgdP = [ diag(dsdx) * dgdP; dsdp' * dgdP];
