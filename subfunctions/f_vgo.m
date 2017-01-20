function [fx] = f_vgo(x,P,u,in)
fx = x + sig(P(1))*(P(2)*u-x);