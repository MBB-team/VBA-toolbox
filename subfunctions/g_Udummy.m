function [gx,dgdx,dgdp] = g_Udummy(x,P,u,in)
gx = P(1) + u(1).*P(2) +u(2).*P(3);
dgdx = [];
dgdp = [1;u(1);u(2)];