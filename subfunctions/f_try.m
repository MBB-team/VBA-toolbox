function [fx,dfdp] = f_try(X,P)
% dummy sigmoidal evolution function
[SP,dsdp] = sigm(P,struct('G0',2,'S0',-1,'beta',1,'INV',0),[0;-1]);
Je = reshape(SP,3,3);
fx = Je*X;
In = eye(3);
dfdp = diag(dsdp)*kron(X,In);