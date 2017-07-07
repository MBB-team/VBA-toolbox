function [gx] = g_sig(x,P,u,in)
gx = 1./(1+exp(-P(1)*x+P(2)));