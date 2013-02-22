function [OBF,ep] = PolyaOBF(n,K,e)
% asymptotic OBF with equal counts
OBF = K.*log((K-1)./(K+n)) + n.*log(n./(K+n)) - log(K-1) +1;

gr = 0:1e-3:1;
a = n/2+1+[e;-e];
pdf = (gr.^(a(1)-1)).*((1-gr).^(a(2)-1));
ind = min(find(gr>=0.5));
pdf = pdf./sum(pdf);
ep = sum(pdf(ind:end));
% cdf = trapz(gr(ind:end),pdf(ind:end));
% lep = gammaln(sum(a)) - sum(gammaln(a)) + log(cdf);


