function P = getFamily(L,families)
[K,n] = size(L);
nf = length(families);
C = zeros(K,nf);
for i=1:nf
    indf = families{i};
    C(indf,i) = 1;
end
f0 = C*sum(C,1)'.^-1/size(C,2);
P = zeros(nf,n);
for i=1:n
    logPm = L(:,i) + log(f0);
    Pm = exp(logPm-max(logPm));
    Pm = Pm./sum(Pm);
    P(:,i) = C'*Pm;
end