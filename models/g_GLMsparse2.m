function [gx,dgdx,dgdP] = g_GLMsparse2(x,P,u,in)
[sP,dsdx,dsdP] = sparsify(P(1:end-1),P(end));
[gx,dgdx,dgdP] = g_GLM(x,sP,u,in);
dgdx = [];
% dgdP = [];
dgdP = [dsdx*dgdP;dsdP'*dgdP];
