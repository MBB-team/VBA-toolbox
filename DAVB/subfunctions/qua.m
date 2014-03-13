function [gx] = qua(x,P,u,in)

gx = P(1) + P(2)*u + P(3)*u.^2;
[gx] = augmentgx(gx,in);
gx = gx(:);