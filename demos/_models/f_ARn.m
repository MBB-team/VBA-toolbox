function [fx,dfdx,dfdp] = f_ARn(x,theta,u,in)

% AR(1) evolution function with exponential decay
n = size(x,1);
[alpha,dsda] = VBA_sigmoid(theta(1:n));
xf = theta(n+1:2*n);
fx = x + -alpha(:).*(x-xf);

dfdx = eye(n) - diag(alpha);
dfdp(1:n,1:n) = -diag((x-xf).*dsda(:));
dfdp(n+1:2*n,1:n) = diag(alpha);