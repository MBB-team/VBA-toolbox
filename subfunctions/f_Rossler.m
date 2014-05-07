function [fx] = f_Rossler(x,P,u,in)

% Rossler oscillator evolution function

a = P(1);
b = P(2);
c = P(3);

xdot = - x(2) - x(3);
ydot = x(1) + a.*x(2);
zdot = b + x(3).*(x(1)-c);

fx = x + in.deltat.*[xdot;ydot;zdot];