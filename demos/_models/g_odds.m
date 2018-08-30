function [gx,dgdx] = g_odds(x,P,u,in)

x = x-max(x);
gx = (exp(x)+0)./sum(exp(x)+0);

n = size(x,1);
I = eye(n);
dgdx = (I - gx*ones(1,n))*diag(gx);
