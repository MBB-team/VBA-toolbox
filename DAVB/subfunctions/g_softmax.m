function  [ gx,dgdx,dgdP ] = g_softmax( x_t,P,u_t,in )

%%% Fonction d'apprentissage
% INPUT
% - x_t : Les états cachés sont les Qvalues (2*1)
% - P : beta du softmax (1*1)
% - u_t : l'action et la récompense ( 2*1 = 1:actions, 2:récompenses)
% - in : []
% OUTPUT
% - gx : P(a=a1|x_t)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

beta = exp(P);
a = 1;
na = 2;

dQ = (x_t(a)-x_t(na));
gx =sig( beta*dQ );

dgdx = zeros(2,1);
dgdx(a) = beta*gx*(1-gx);
dgdx(na) = -beta*gx*(1-gx);

dgdP = [beta*dQ*gx*(1-gx)];


function y=sig(x)
y = 1/(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;