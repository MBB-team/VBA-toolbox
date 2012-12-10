function x = Av(vecA,v)

n = sqrt(length(vecA));
A = reshape(vecA,n,n);
x = A*v;