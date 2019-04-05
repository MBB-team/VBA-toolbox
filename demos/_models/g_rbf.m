function [gx] = g_rbf(x,P,u,in)

centers = in.centers + exp(P(1));
sig = exp(P(2));

N = length(in.centers);
Xrbf = zeros(length(in.grid),N);
for i=1:N
    Xrbf(:,i) = exp(-0.5*(centers(i)-in.grid).^2./sig);
    Xrbf(:,i) = Xrbf(:,i)./sum(Xrbf(:,i));
end

gx = zeros(N,1);
for i=1:N
    y = Xrbf(:,i);
    X = Xrbf(:,setdiff(1:N,i));
    P0 = eye(length(in.grid)) - X*pinv(X'*X)*X';
    err = y'*P0*y;
    gx(i) = sum(err);
end

gx = gx + P(3);

if ~in.corr
    gx = Xrbf;
end