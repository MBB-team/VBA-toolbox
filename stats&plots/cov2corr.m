function [C] = cov2corr(V)
% computes correlation matrix from covariance matrix

dv = diag(V);
n = size(V,1);
C = zeros(n,n);
for i=1:n
    for j=1:n
        C(i,j) = V(i,j)./sqrt(dv(i).*dv(j));
    end
end

