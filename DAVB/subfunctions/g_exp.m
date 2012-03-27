function [gx] = g_exp(x,phi,u,in)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
gx = phi(2)*exp(phi(1)*in.x)+phi(3);