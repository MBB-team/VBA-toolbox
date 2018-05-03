function [fx,dF_dX,dF_dTheta] = f_Lorenz(Xt,Theta,ut,inF)
% Lorenz chaotic evolution function

deltat = inF.deltat;

x       = Xt;
rho     = Theta(1);
sigma   = Theta(2);
beta    = Theta(3);

fx      = zeros(3,1);
fx(1)   = sigma.*(x(2)-x(1));
fx(2)   = x(1).*(rho-x(3))-x(2);
fx(3)   = x(1).*x(2) - beta.*x(3);
fx      = deltat.*fx + x;

J       = zeros(3,3);
J(1,:)  = [-sigma,sigma,0];
J(2,:)  = [rho-x(3),-1,-x(1)];
J(3,:)  = [x(2),x(1),-beta];
dF_dX   = deltat.*J' + eye(3);

dF_dTheta = zeros(3,3);
dF_dTheta(1,:)      = [0,x(1),0];
dF_dTheta(2,:)      = [x(2)-x(1),0,0];
dF_dTheta(3,:)      = [0,0,-x(3)];
dF_dTheta           = deltat.*dF_dTheta;
                   