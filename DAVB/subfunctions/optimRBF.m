function [Xrbf,in] = optimRBF(T,N,dt)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
t = 0:dt:T;
rbfdt = T./N;
in.centres = [0:N-1]*rbfdt + 0.5*rbfdt;
in.sig = 0.25*rbfdt;

[Xrbf,dudx,dudp] = u_RBF([],eye(N),t,in);



