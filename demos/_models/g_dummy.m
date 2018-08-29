function [gx] = g_dummy(x,P,u,in)
gx = 1./(1+exp(-x./P));
