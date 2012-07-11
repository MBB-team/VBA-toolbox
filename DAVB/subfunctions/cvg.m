function [gx] = cvg(x,P,u,in)

gx = P(1) + P(2)*(1-exp(-exp(P(3))*u));
[gx] = augmentgx(gx,in);
gx = gx(:);