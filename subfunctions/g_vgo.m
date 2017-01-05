function [gx] = g_vgo(x,P,u,in)
gx = sig(exp(P(1))*x+P(2));