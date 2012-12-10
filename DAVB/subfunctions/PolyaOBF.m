function OBF = PolyaOBF(a0,nk)
K = size(nk,1);
n = sum(nk);
lOBF = -gammaln(a0+n);
for k=1:K
    lOBF = lOBF+gammaln(a0+nk(k))+nk(k)*log(K);
end
OBF = exp(lOBF);

