function  [ gx ] = g_probability_matching( x_t,P,u_t,in )
%%% Probability of emission for TWO alternatives with CONTINUOUS outcome
% There are two options (1 and 2)
% For each option there is a belief about the associated outcome.
% This has the form of a Gaussian distribution 
% Subject has a utility function u.
% Decision is made through sampling from the probability of the
% decision to lead to the highest reward (hence the name: probability
% matching)

% INPUT
% - x_t : 4*1 vector : sufficient statistics of outcome for each alternative [mu1;s1;mu2;s2]
% - P : []
% - u_t : []
% - in.u : handle of the utility function
% OUTPUT
% - gx : Scalar, P(a=a1|x_t), probability of action 1

mu1 = x_t(1);
mu2 = x_t(3);
s1 = sqrt(x_t(2));
s2 = sqrt(x_t(4));
[du1,u1] = numericDiff(in.u_fname,1,mu1,P,u,in);
[du2,u2] = numericDiff(in.u_fname,1,mu2,P,u,in);
% [u1,du1]=feval(in.u_fname,mu1); % utility at mean %numericdiff
% [u2,du2]=feval(in.u_fname,mu2); % utility at mean

sU1 = sqrt(du1^2*s1^2);
sU2 = sqrt(du2^2*s2^2);
gx  = sig((u2-u1)/(sqrt(3)*sU1/pi)/sqrt(1+sU2^2/sU1^2));

function y=sig(x)
y = 1/(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps; 
