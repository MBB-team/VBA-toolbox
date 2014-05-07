function [fx,dF_dX,dF_dTheta] = f_Henon(Xt,Theta,ut,inF)
% Henon chaotic map evolution function

x       = Xt;
a       = Theta(1);
b       = Theta(2);

fx      = zeros(2,1);
fx(1)   = x(2) + 1 -a.*x(1).^2;
fx(2)   = b.*x(1);

J       = zeros(2,2);
J(1,:)  = [-2.*a*x(1),1];
J(2,:)  = [b,0];
dF_dX   = J';

dF_dTheta = zeros(2,2);
dF_dTheta(1,:)      = [-x(1).^2,0];
dF_dTheta(2,:)      = [0,x(1)];
dF_dTheta = dF_dTheta';


                       