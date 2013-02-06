function  u = visualU(P,in)


A1 = in.A1;
A2 = in.A2;
mu3_0 = in.mu3_0;
[n1,n2] = size(A1);
n3 = size(A2,2);
s0 = P(1);
s1 = P(2);
s2 = P(3);
s3 = P(4);
x3 = ones(size(mu3_0));
x2 = A2*x3 + 0*s2*randn(n2,1);
x1 = A1*x2 + 0*s1*randn(n1,1);
u = x1 + s0*randn(n1,1);
u = u.*repmat([1;0],n1/2,1);