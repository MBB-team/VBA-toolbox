function [gx,dgdx,dgdP] = g_GLMsparse(x,P,u,in)
[sP,dsdP] = sparseTransform(P,in.sparseP);
[gx,dgdx,dgdP] = g_GLM(x,sP,u,in);
dgdP = dsdP*dgdP; % for exploiting the analytical gradients from g_GLM
