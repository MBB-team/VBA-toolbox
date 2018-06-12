function [fx] = f_2d(x,P,u,in)
xdot = in.A*x;
fx = x + in.dt*xdot;
