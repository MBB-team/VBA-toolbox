function [gx,dgdx] = g_RFX(x,P,u,in)

gx = in.X*x;
dgdx = in.X';