function [g] = dsdv(x,P)
try % parameters of the Gaussian observation filter
    P.sig;
catch
    P.sig = struct('r',0.54,'eta',0);
end
sx = 1./(1+exp(P.sig.r.*(P.sig.eta-x)));
g = P.sig.r.*sx.*(1-sx);