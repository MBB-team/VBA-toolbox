function [gx,dgdx,dgdp] = g_classif(x,P,u,in)

gx = sss([in.X',ones(size(in.X,2),1)]*P);
dgdx = [];
dgdp = diag(gx.*(1-gx))*[in.X',ones(size(in.X,2),1)];
dgdp = dgdp';

function sx = sss(x)
sx = 1./(1+exp(-x));
sx(sx < 1e-8) = 1e-8;
sx(sx > 1-1e-8) = 1-1e-8;

