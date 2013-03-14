function [x] = mexicanhat(lx,rx,n,sigma)

t2 = (lx:(rx-lx)/(n-1):rx).^2;

x = (2/(sqrt(3*sigma)*pi^.25)) * (1-t2/sigma^2).*exp(-t2/(2*sigma^2));

end