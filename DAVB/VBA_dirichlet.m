function p = VBA_dirichlet(a,x)

lp = sum((a-1).*log(x)) + gammaln(sum(a)) - sum(gammaln(a));
p = exp(lp);
