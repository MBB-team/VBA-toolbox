function [gx] = pow(x,P,u,in)

gx = P(1) + P(2).*u.^P(3);
[gx] = augmentgx(gx,in);
gx = gx(:);