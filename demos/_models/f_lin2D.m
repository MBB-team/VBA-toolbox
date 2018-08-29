function [fx,dF_dX,dF_dTheta] = f_lin2D(Xt,Theta,ut,inF)
% dummy 2D linear evolution function

deltat = inF.deltat;

if ~ isfield (inF, 'a')
    inF.a = 1;
end
if ~ isfield (inF, 'b')
    inF.b = 1e-1;
end

inF.a = inF.a * exp (Theta(1));

A = [- inF.b, - inF.a;
     1, - inF.b];
    
fx = Xt + deltat * (A * Xt + ut);
dF_dX = eye (size (Xt, 1)) + deltat * A';

dF_dTheta = deltat * [- inF.a * Xt(2), 0];