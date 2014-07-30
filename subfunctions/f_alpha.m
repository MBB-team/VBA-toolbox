function [fx] = f_alpha(x,P,u,in)
H = P(1);
T = P(2);
fx = x + in.dt.* [x(2);(H/T)*u-(2/T)*x(2)-(1/T^2)*x(1)];