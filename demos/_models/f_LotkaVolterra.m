function [fx] = f_LotkaVolterra(X,P,u,in)

A = reshape(P,3,3);

xdot = diag(X)*A*(X-ones(3,1));
fx = X + in.deltat*xdot;
