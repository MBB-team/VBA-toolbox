function [gx,dgdx] = g_RFX(x,P,u,in)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

gx = in.X*x;
dgdx = in.X';