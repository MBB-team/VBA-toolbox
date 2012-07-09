function  [ gx ] = g_softmax_utility( x_t,P,u_t,in )
%%% Softmax emission for TWO alternatives with CONTINUOUS outcome
% There are two options (1 and 2)
% For each option there is a belief about the associated outcome.
% This has the form of a Gaussian distribution 
% Subject has a utility function u.
% Decision is made by softmaximizing the expected utility of both actions

% INPUT
% - x_t : 4*1 vector : sufficient statistics of outcome for each alternative [mu1;s1;mu2;s2]
% - P : Scalar, inverse temperature of the softmax decision rule
% - u_t : []
% - in.u : handle of the utility function
% OUTPUT
% - gx : Scalar, P(a=a1|x_t), probability of action 1

mu1 = x_t(1);
mu2 = x_t(3);
[u1,du1]=feval(in.u_fname,mu1) % utility at mean
[u2,du2]=feval(in.u_fname,mu2) % utility at mean
dQ = (u1-u2);
b= exp(P(1));
gx =sig( -b*dQ );

function y=sig(x)
y = 1/(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps; 
