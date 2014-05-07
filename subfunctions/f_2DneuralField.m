function [fx,dF_dX,dF_dTheta] = f_2DneuralField(Xt,Theta,ut,inF)
% 2D neural field evolution function (wave propagation)
% function [fx,dF_dX,dF_dTheta,d2F_dXdTheta] = f_2DneuralField(Xt,Theta,ut,inF)
% This function evaluates the evolution function of a linear 2D neural
% field.

deltat = inF.deltat;
try
    L = inF.L;
catch
    n = sqrt(size(Xt,1)/2);
    I = speye(n,n);
    E = sparse(2:n,1:n-1,1,n,n);
    D = E+E'-2*I;
    L = kron(D,I)+kron(I,D);
end

K = Theta(1);
c = Theta(2);
b = Theta(3);


n = size(Xt,1);
ns = sqrt(n./2);
x1 = Xt(1:n/2);
x2 = Xt(n/2+1:end);

Lx = L*x1;

f1 = zeros(n/2,1);
f1 = x2;
f2 = zeros(n/2);
f2 = -2*K*x2 - K.^2.*x1 + 1.5*c^2.*Lx + b.*ut;

fx = Xt + deltat.*[f1;f2];

In = speye(n/2);
On = zeros(n/2,1);

dfdx = [    zeros(n/2)              In
           -K.^2.*In + c^2.*L      -2*K*In  ]; 
dF_dX = eye(n) + deltat.*dfdx';

dfdp = zeros(n,3);
dfdp(:,1) = [ On ; -2*x2 - 2*K.*x1];
dfdp(:,2) = [ On ; 2*c.*Lx];
dfdp(:,3) = [ On ; ut];

dF_dTheta = deltat.*dfdp';




