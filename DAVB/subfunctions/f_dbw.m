function [fx] = f_dbw(x,theta,u,in)

% dummy double-well evolution function
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

dt = in.dt;

xdot = ((2*(4*x - 8))/exp(2*(x - 2)^2) + (x/4 + 1/2)/exp((x + 2)^2/8))/(2/exp(2*(x - 2)^2) + 1/exp((x + 2)^2/8));
fx = x - dt*xdot;