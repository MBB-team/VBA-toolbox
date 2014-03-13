function [gx] = cvg2(x,P,u,in)

gx = P(1) + P(2)*exp(-exp(P(3))*u);
[gx] = augmentgx(gx,in);
gx = gx(:);