function gx = cst(x,P,u,in)

gx = P(1).*ones(size(u));
[gx] = augmentgx(gx,in);
gx = gx(:);