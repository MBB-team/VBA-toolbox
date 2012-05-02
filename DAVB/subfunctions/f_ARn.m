function [fx] = f_ARn(x,theta,u,in)

% AR(1) evolution function with exponential decay
n = size(x,1);
in.G0 = 1;
in.S0 = 0;
in.beta = 1;
in.INV = 0;
[alpha,dsda] = sigm(theta(1:n),in,[]);
xf = theta(n+1:2*n);
fx = x + -alpha(:).*(x-xf);
% 
% dfdx = eye(n) - diag(alpha);
% dfdp(1:n,1:n) = -diag((x-xf).*alpha);
% dfdp(n+1:2*n,1:n) = diag(alpha);