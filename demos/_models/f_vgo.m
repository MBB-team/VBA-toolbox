function [fx] = f_vgo(x,P,u,in)
fx = x + VBA_sigmoid(P(1))*(P(2)*u-x);