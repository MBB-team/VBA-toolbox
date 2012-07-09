function [gx] = g_Mathys(x,P,u,in)
% computes the probability of the learner choosing 'a=1'.
% [gx] = g_VBvolatile0(x,P,u,in)
% This function derives the action emission law from the current belief of
% the learner. It is valid in an empirical context akin to a multiarmed
% bandit problem, whereby the learner is asked to bet on the outcome at
% each trial.
% Note: the action emission law accounts for uncertainty in the hidden
% states x2, when passed through the sigmoid mapping.
% IN:
%   - x: the previous posterior sufficient statistics.
%   - P: the response model parameters vector.
%   - u: the current input to the learner.
%   - in: further quantities handed to the function.
% OUT:
%   - gx: the expected next input to the learner, according to its previous
%   updated belief about the hidden states.


p1 = sgm(x(2),1);
gx = p1;

%gx = sgm(exp(P(1)*(1-2*p1),1)); % Softmax decision on probability

% P(1): temperature/utility scaling
% P(2): preference bias towards choosing action a=1