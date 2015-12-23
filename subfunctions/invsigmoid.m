function   y = invsigmoid(x) % pour x entre 0 et 1)

y=log(x./(1-x));
ind0=find(x<=1e-14);
ind1=find(x>=(1-1e-14));
y(ind0)=-10^4;
y(ind1)=10^4;
end
