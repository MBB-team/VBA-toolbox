function [F,F0,a] = PolyaLEV(m,a0)
[K,n] = size(m);
a = m*ones(n,1) + a0;
F = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(a)) - gammaln(n+K*a0);
F0 = -n*log(K);