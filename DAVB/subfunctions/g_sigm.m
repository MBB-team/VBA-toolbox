function [Sx,dsdx,dsdp] = g_sigm(x,Phi,u,in)
% sigmoid observation mapping
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
try; in.x; catch; in.x = 1; end
if in.x
    [Sx,dsdx,dsdp] = sigm(x,in,Phi);
else
    [Sx,dsdp,dsdx] = sigm(Phi,in,x);
    dsdp = diag(dsdp);
end