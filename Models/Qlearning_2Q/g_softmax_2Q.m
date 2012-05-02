function  [ gx ] = g_softmax_2Q( x_t,P,u_t,in )
%%% Learning function for 2 Qvalues
% INPUT
% - x_t : 4*1 vector, [Qpost1;Qpre1;Qpost2;Qpre2]
% - P : Scalar, inverse temperature of the softmax decision rule
% - u_t : 2*1 vector, Action and reward
% OUTPUT
% - gx : Scalar, P(a=a1|x_t), probability of action 1
%gx = feval(...
%   @g_softmax,... % the function to be called
%   x_t([2,4]),... % the indices of the hidden states concerned by the evolution function 
%   P(1),... % the indices of the evolution parameters
%   u_t,... % the indices of the data
%   []); 
    
beta = exp(P);
dQ = (x_t(1)-x_t(2));
gx =sig( -beta*dQ );


function y=sig(x)
y = 1/(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps; 



% Computing derivatives
% dgdx = zeros(2,1);
% dgdx(a) = beta*gx*(1-gx);
% dgdx(na) = -beta*gx*(1-gx);
% 
% dgdP = [beta*dQ*gx*(1-gx)];ontributions are added