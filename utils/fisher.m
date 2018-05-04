function f = fisher(r,inv)
% apply Fisher transformation to Pearson correlation coefficient
try,inv;catch;inv=0;end
if ~inv
    f = 0.5*log((1+r)./(1-r));
else
    f = (exp(2*r)-1)./(1+exp(2*r));
end