function [GX,dGdX,dGdP] = g_embedAR(X,P,ut,in)
% AR(1) evolution function

n = size(X,1);
x = X(1:n/2);
z = X(n/2+1:n);

[gx,J,dgdP] = VBA_evalFun('g',x,P,ut,in.opt,in.dim,0);

GX = gx;
dGdX = [J;zeros(n/2)];
dGdP = dgdP;
