function [gx] = g_nl0(x,P,u,in)

try
    gx = P(1)*exp(P(2)*in.x);
catch
    gx = x(1)*exp(x(2)*P);
end