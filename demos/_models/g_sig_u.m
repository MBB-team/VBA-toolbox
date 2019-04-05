function [Sx] = g_sig_u(x,Phi,u,in)
% sigmoid observation mapping
[Sx] = VBA_sigmoid(u,'slope',exp(Phi(1)),'center',Phi(2));