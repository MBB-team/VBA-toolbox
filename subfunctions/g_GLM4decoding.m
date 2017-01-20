function [gx,dgdx,dgdp] = g_GLM4decoding(Xt,P,ut,in)


Xt = Xt+P(1);
dgdx = 1;
dgdp = 1;
