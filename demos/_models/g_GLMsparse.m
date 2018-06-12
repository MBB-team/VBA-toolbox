function [gx,dgdx,dgdP] = g_GLMsparse(x,P,u,in)
[sP, dsdP] = VBA_sparsifyPrior (P);
[gx,dgdx,dgdP] = g_GLM(x,sP,u,in);
dgdP = diag(dsdP)*dgdP; % for exploiting the analytical gradients from g_GLM
