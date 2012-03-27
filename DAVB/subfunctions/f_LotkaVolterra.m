function [fx] = f_LotkaVolterra(X,P,u,in)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

A = reshape(P,3,3);

xdot = diag(X)*A*(X-ones(3,1));
fx = X + in.deltat*xdot;
