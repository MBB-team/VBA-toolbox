function  [ gx,dgdx,dgdP ] = g_softmax( x_t,P,u_t,in )
% INPUT
% - x_t : Q-values (2*1)
% - P : inverse temperature (1*1)
% - u_t : previous action and feedback
% - in : []
% OUTPUT
% - gx : P(a=a1|x_t)

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