function [I,NextSigma,NextdeltaMu] = VarVolatility(x,mu,va,a,b)
% OTO: varational energy (and curvature) of volatile observer

if ~isinf(va)
    iva = va.^-1;
else
    iva = 0;
end
ex = exp(x(:));
v = a + ex;
I = -0.5.*iva*(x(:)-mu).^2 ...
    -0.5*log(v) ...
    -0.5*b./v;
dIdx = iva*(mu-x(:)) - 0.5.*ex./v + 0.5.*b.*ex./v.^2;
dI2dx2 = -iva -0.5.*a*ex./v.^2 -0.5.*b.*ex.*(ex-a)./v.^3;
NextSigma = -dI2dx2.^-1;
NextdeltaMu = NextSigma.*dIdx;




