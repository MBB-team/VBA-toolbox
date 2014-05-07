function [FX,dFdX,dFdP] = f_embedAR(X,P,ut,in)
% AR(1) embedding evolution function

n = size(X,1);
x = X(1:n/2);
z = X(n/2+1:n);
% call native evolution function
[fx,J,dfdP] = VBA_evalFun('f',x,P,ut,in.opt,in.dim,0);
% add AR(1) noise and form derivatives
FX = [fx+z;z];
dFdX = [ J          zeros(n/2)
         eye(n/2)   eye(n/2) ];
dFdP = [dfdP,zeros(size(P,1),n/2)];
