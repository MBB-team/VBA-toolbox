function [fx,dfdx,dfdp] = f_gen(Xt,Theta,ut,inF)

% Generic evolution function (up to quadratic terms)

% [fx,dfdx,dfdp] = f_gen(Xt,Theta,ut,inF)

deltat = inF.deltat;
n = size(Xt,1);
nc = factorial(n)./(factorial(2)*factorial(n-2));
A = reshape(Theta(1:n^2),n,n);
B = reshape(Theta(n^2+1:n*(2*n+nc)),n,n+nc);

xij = zeros(n+nc,1);
ind = cell(n,1);
ind2 = zeros(n,1);
k = 0;
for i=1:n
    for j=1:n
        if j >= i
            k = k+1;
            xij(k) = Xt(i).*Xt(j);
            if i == j
                ind2(i) = k;
            else
                ind{i} = [ind{i},k];
                ind{j} = [ind{j},k];
            end
        end
    end
end
dbxdx = zeros(n,n);
for i=1:n
    for j=1:n
        dbxdx(i,j) = 2*B(i,ind2(j)).*Xt(j) ...
            + B(i,ind{j})*Xt(setdiff(1:n,j));
    end
end

f = A*Xt + B*xij;
fx = Xt + deltat.*f;
dfdx = eye(n) + deltat*(A + dbxdx)';

dfdp = zeros(n,n*(2*n+nc));
dfdp(:,1:n^2) = kron(Xt',eye(n));
dfdp(:,n^2+1:n*(2*n+nc)) = kron(xij',eye(n));
dfdp = deltat*dfdp';

