function y = sgm(x,a)
% simple sigmoid mapping
% function y = sgm(x,a)
if nargin == 1
    a = 1;
end
localEps = 1e-4;

y = a./(1+exp(-x));
y(y<a*localEps)=a*localEps;
y(y>a.*(1-localEps))=a.*(1-localEps);
