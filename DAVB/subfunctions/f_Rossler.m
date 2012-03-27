function [fx] = f_Rossler(x,P,u,in)
% Rossler oscillator evolution function
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

a = P(1);
b = P(2);
c = P(3);

xdot = - x(2) - x(3);
ydot = x(1) + a.*x(2);
zdot = b + x(3).*(x(1)-c);

fx = x + in.deltat.*[xdot;ydot;zdot];

% J = zeros(3,3);
% J(1,:) = [0,-1,-1];
% J(2,:) = [1,a,0];
% J(1,:) = [x(3),0,x(1)-c];
% dfdx = eye(3) + in.deltat.*J';
% 
% dfdp = zeros(3,3);
% dfdp(1,:) = zeros(1,3);
% dfdp(2,:) = [x(2),0,0];
% dfdp(3,:) = [0,1,-x(3)];
% dfdp = in.deltat.*dfdp';