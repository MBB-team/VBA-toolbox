function [Ux,dUdx,dUdP,d2Udx2,d2UdP2] = U_dummy(x,P,u,y,in)
Ux = -0.5*(y-u(1).*P(2)-u(2).*P(3)-P(1)).^2;
dUdx = [];
dUdP = (y-u(1).*P(2)-u(2).*P(3)-P(1))*[1;u(1);u(2)];
d2Udx2 = [];
d2UdP2 = -[ 1       u(1)    u(2)
            u(1)    u(1).^2 u(1)*u(2)
            u(2)    u(1)*u(2)    u(2).^2 ];