function [fx,dfdx] = f_ARplus(x,P,u,in)
% AR(1) evolution function
I = eye(length(x));
fx = (I+diag(P))*x;
dfdx = I+diag(P);