function [g] = dsdv(x,P)
try % parameters of the Gaussian observation filter
    P.sig;
catch
    P.sig = struct('r',0.54,'eta',0,'g',0.135);
end
g = P.sig.r.*exp(P.sig.r.*(P.sig.eta-x))./(1+exp(P.sig.r.*(P.sig.eta-x))).^2;