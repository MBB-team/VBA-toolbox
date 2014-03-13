function [G] = dispersivePropagator(r,t,P)

H = zeros(size(t));
H(t>0) = 1;

v = P(1);
s = P(2);
G0 = P(3);

if r~=0
    G = (t.^-1).*exp(-(r.^2+(v*t).^2)./(2*s*v*t)).*H;
else
    G = (t.^-1).*exp(-(v*t)./(2*s)).*H;
end
G = (G0./(4*pi*s^2))*G;