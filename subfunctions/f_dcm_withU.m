function [fx,dF_dX,dF_dP] = f_dcm_withU(x,P,ut,in)
% linear evolution function with parameterized input

deltat = in.deltat;

n = size(x,1);
A = reshape(P(in.indA),n,n);

if ~isempty(in.u_fname)
    [uu,dudx,dudp] = feval(in.u_fname,[],P(in.indU),ut(n+1),in);
else
    uu = ut(1:n);
    dudp = zeros(length(in.indU),n);
end

fx = x + deltat.*(A*x + uu(:));
dF_dX = eye(n) + deltat.*A';

dfdp = kron(x,eye(n));
dfdp = [dfdp;dudp];
dF_dP = deltat*dfdp;