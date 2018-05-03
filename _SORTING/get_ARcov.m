function Sigma = get_ARcov(n)
I = eye(n);
D = diag(ones(n-1,1),-1);
Sigma = (I-D)*(I-D)';
Sigma(n,n) = 1;