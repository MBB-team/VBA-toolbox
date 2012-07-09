function [gx] = g_softmax_EU_2p(x,P,u,in)
% computes the probability of the learner choosing 'a=1' through
% a softmax decision based on the expected utility

% [gx] = g_VBvolatile1(x,P,u,in)
% This function derives the action emission law from the current belief of
% the learner. It is valid in an empirical context akin to a multiarmed
% bandit problem, whereby the learner is asked to bet on the outcome at
% each trial.

% IN:
%   - x: the previous posterior sufficient statistics.
%   - P: the response model parameters vector.
%   - u: the current input to the learner.
%   - in: further quantities handed to the function.
%       - u_fname: handle of the utility function
% OUT:
%   - gx: the expected next input to the learner, according to its previous
%   updated belief about the hidden states.

p1_1 = sgm(x(2),1); % probability of outcome 1 given choice 1
p1_2 = sgm(x(7),1); % probability of outcome 1 given choice 2
    
[u1]=feval(in.u_fname,in.o1,P(2),in); % utility of outcome 1
[u2]=feval(in.u_fname,in.o2,P(2),in); % utility of outcome 2
    
gx  = sig(-(u1-u2)*(p1_1-p1_2)*exp(P(1)));

function y=sig(x)
y = 1/(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps; 
