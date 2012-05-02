function [fx,dF_dX,dF_dTheta] = f_lin2D(Xt,Theta,ut,inF)
% dummy 2D linear evolution function

deltat = inF.deltat;

try
    a = inF.a;
catch
    a = 1;
end
try
    b = inF.b;
catch
    b = 1e-1;
end

a = a.*exp(Theta(1));

A = [   -b  -a
        1   -b];
    
fx = Xt + deltat.*A*Xt;
dF_dX = eye(size(Xt,1)) + deltat.*A';

dF_dTheta = deltat.*[-a.*Xt(2),0];