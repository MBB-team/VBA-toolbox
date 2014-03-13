function [I] = optDT(x,dF,Phi,fname)
% generic posterior risk for reaction times
alpha = Phi(1);
cdot = feval(fname,x,Phi(2:end));
I = -(exp(-alpha.*x)-(cdot).^2./(alpha.*dF)).^2;
