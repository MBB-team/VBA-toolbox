function  [ gx,dgdx,dgdP ] = g_softmax(x,P,u,in )
% INPUT
% - x : Q-values (2x1)
% - P : inverse temperature (1x1)
% - u : [useless]
% - in : [useless]
% OUTPUT
% - gx : P(a=1|x)

beta = exp(P);
dQ = (x(1)-x(2));
gx = sig( beta*dQ );
dgdx = zeros(size(x,1),1);
dgdx(1) = beta*gx*(1-gx);
dgdx(2) = -beta*gx*(1-gx);
dgdP = [beta*dQ*gx*(1-gx)];
