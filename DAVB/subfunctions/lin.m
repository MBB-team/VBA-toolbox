function [gx] = lin(x,P,u,in)

gx = P(1) + P(2)*u;
[gx] = augmentgx(gx,in);
gx = gx(:);