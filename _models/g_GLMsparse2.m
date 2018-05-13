function [gx,dgdx,dgdP] = g_GLMsparse2(x,P,u,in)


[sP, dsdx, dsdP] = VBA_sparsifyPrior (P(1:end-1), P(end));

[gx, ~, dgdP] = g_GLM (x, sP, u, in);

dgdx = [];

dgdP = [ diag(dsdx) * dgdP; dsdP' * dgdP];
