function [gx] = g_vgo(x,P,u,in)
gx = sig(exp(P)*x);