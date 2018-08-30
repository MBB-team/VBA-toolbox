function [fx] = f_2dwu(x,P,u,in)
xdot = in.A*x + diag(P)*u;
fx = x + in.dt*xdot;