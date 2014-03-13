function [fx,dfdx,dfdp] = f_gradf(X,Theta,u,in)
% gradient of a dummy quadratic function
fdot = -exp(Theta(1))*(X-Theta(2));

fx = X + fdot;

dfdx = 1 -exp(Theta(1));
dfdp = [fdot;exp(Theta(1))];