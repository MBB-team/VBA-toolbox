function [fx] = f_vgo(x,P,u,in)
fx = x + sig(P)*(u-x);