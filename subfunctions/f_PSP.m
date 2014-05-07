function [fx,dfdx] = f_PSP(x,P,u,in)
% evolution function for postsynaptic potentials

H = P(1);
Krise = 1/P(2);
if length(P) == 2
    Kdecay = Krise;
else
    Kdecay = 1/P(3);
end

J = [  0                1
       -Krise*Kdecay    -2*Krise  ];

xdot = J*x + [0;1]*H*Krise*u;   
fx = x + in.dt*xdot;
dfdx = eye(2) + in.dt*J;
