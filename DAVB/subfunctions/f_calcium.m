function [fx,dF_dX,dF_dTheta] = f_calcium(Xt,Theta,ut,inF)

% (2-photons) calcium imaging evolution function

deltat = inF.delta_t;
try
    a = inF.a.^-1;
catch
    a = 1;
end
a = a.*exp(Theta(1));
b = exp(Theta(2));
tmp = [-a*Xt(1)+exp(Xt(2));-b*Xt(2)];

fx = Xt + deltat.*tmp + [1;0]*ut;

A = [ -a  exp(Xt(2))
       0  -b          ];
dF_dX = eye(2) + deltat.*A';

dF_dTheta = [ -deltat.*Xt(1)*a  0
              0                 -deltat.*Xt(2)*b];