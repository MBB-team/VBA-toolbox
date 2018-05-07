function [Sx,dsdx,dsdp] = g_sigm(x,Phi,u,in)
% sigmoid observation mapping
try
    in.x; 
catch
    in.x = 1;
end

if in.x
    [Sx,dsdx,dsdp] = VBA_sigmoid(x);
else
    [Sx,dsdp,dsdx] = sigm(Phi,in,x);
    dsdp = diag(dsdp);
    pause()
end