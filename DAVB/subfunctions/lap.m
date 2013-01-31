function L = lap(n)

Un = ones(n-1,1);
L = -diag(Un,-1) - diag(Un,1) + 2*eye(n);
L(1,1) = 1;
L(n,n) = 1;