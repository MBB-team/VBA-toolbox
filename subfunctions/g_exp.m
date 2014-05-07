function [gx] = g_exp(x,phi,u,in)

gx = phi(2)*exp(phi(1)*in.x)+phi(3);