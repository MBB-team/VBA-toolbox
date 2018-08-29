function [fx, dfdx, dfdp] = f_Rossler(x,P,u,in)

% Rossler oscillator evolution function

a = P(1);
b = P(2);
c = P(3);

xdot = - x(2) - x(3);
ydot = x(1) + a.*x(2);
zdot = b + x(3).*(x(1)-c);

fx = x + in.deltat.*[xdot;ydot;zdot];

dfdx = eye(3) + in.deltat * [0, 1, x(3); - 1, a, 0; -1, 0, (x(1)- c)];

dfdp = in.deltat * [0, x(2), 0; 0, 0, 1; 0, 0, - x(3)] ;