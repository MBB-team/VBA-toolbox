function [fx,dfdx,dfdP,d2fdxdP] = g_Id_phi(x,P,u,in)
% this function applies the identity mapping
fx = P+0.5;
dfdx = [];
dfdP = eye(size(P,1));
d2fdxdP = [];

