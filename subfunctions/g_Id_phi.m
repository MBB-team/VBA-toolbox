function [fx,dfdx,dfdP] = g_Id_phi(x,P,u,in)
% this function applies the identity mapping
fx = P;
dfdx = [];
dfdP = eye(size(P,1));

