function [g,dgdx,dgdp] = g_matmap(x,P,u,in)

nu = size(u,1);
imat = in.dim.p*nu;
matP = reshape(P(1:imat),in.dim.p,nu);
g = P(imat+1:imat+in.dim.p) + matP*u;
dgdx = [];
dgdp = [kron(u,eye(in.dim.p));eye(in.dim.p)];
