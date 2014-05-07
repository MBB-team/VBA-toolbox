function [fx,dfdx,dfdP] = f_Id(x,P,u,in)
% this function applies the identity mapping
fx = x;
dfdx = eye(size(x,1));
dfdP = zeros(size(P,1),size(x,1));