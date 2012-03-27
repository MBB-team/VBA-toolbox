function [Sx,dsdx,dsdp] = g_sigm_binomial(x,Phi,u,in)
% evaluates the sigmoid function for binomial data analysis
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

in.G0 = 1;
in.S0 = 0;
in.beta = 1;
in.INV = 0;
try,in.x;catch,in.x = 0;end
if in.x % for learning effects (sigmoid parameters evolve over time)
    [Sx,dsdx,dsdp] = sigm(u,in,x);
    dsdx = dsdp;
    dsdp = [];
else
    [Sx,dsdx,dsdp] = sigm(u,in,Phi);
    dsdx = [];
end


