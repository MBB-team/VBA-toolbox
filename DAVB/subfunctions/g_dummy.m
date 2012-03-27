function [gx] = g_dummy(x,P,u,in)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
gx = in.X*P(1);
try
    gx = gx + in.X*P(2);
end
