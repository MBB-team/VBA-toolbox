function [fx] = f_dbw(x,theta,u,in)

% dummy double-well evolution function

dt = in.dt;

xdot = ((2*(4*x - 8))/exp(2*(x - 2)^2) + (x/4 + 1/2)/exp((x + 2)^2/8))/(2/exp(2*(x - 2)^2) + 1/exp((x + 2)^2/8));
fx = x - dt*xdot;