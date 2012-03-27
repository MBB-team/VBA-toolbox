function  [ gx ] = g_softmax( x_t,P,u_t,in )
% ,dgdx,dgdP
%%% Fonction d'apprentissage
% INPUT
% - x_t : Les états cachés sont les Qvalues (2*1)
% - P : beta du softmax (1*1)
% - u_t : l'action et la récompense ( 2*1 = 1:actions, 2:récompenses)
% - in : []
% OUTPUT
% - gx : P(a=a1|x_t)

beta = exp(P(1));
gx =sig( -beta*(x_t(1)-x_t(2)) );

% Computing derivatives
% dgdx = zeros(2,1);
% dgdx(a) = beta*gx*(1-gx);
% dgdx(na) = -beta*gx*(1-gx);
% 
% dgdP = [beta*dQ*gx*(1-gx)];


function y=sig(x)
y = 1/(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;