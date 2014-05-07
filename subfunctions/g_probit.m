function [gx] = g_probit(x,P,u,in)
gx = 1./(1+exp(-exp(P)*(x-in.bias)));