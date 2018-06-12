function [gx] = g_odds2(x,P,u,in)

% P = P-max(P);
gx = (exp(P)+0)./sum(exp(P)+0);
% dgdx = [];
% n = size(P,1);
% I = eye(n);
% dgdP = (I - gx*ones(1,n))*diag(gx);