function [Sx] = g_sig_u(x,Phi,u,in)
% sigmoid observation mapping
[Sx] = sigm(u,struct('mat',1),Phi);
