function [F,F0,a] = PolyaLEV(m,a0)
[K,n] = size(m);
a = m*ones(n,1) + a0;
sa = sum(a);
F = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(a)) - gammaln(sa);
F0 = -log(K)*sum(m(:));