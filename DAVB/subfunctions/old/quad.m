function [fx,J] = quad(Xt,theta)
% emulates quadratic terms in nonlinear DCMs
n = size(Xt,1);
In = speye(n);
dxD = zeros(n,n);
dxD2 = dxD;
for i=1:n
    D{i} = reshape(theta((i-1)*n.^2+1:i*n.^2),n,n);
    tmp = Xt(i)*D{i};
    dxD = dxD + tmp;
    tmp(:,i) = tmp(:,i)+D{i}*Xt;
    dxD2 = dxD2 +tmp;
end

fx = dxD*Xt;
J = dxD2';