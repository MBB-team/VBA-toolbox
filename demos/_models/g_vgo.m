function [gx] = g_vgo(x,P,u,in)
gx = VBA_sigmoid(exp(P(1))*x+P(2));