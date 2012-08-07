function [gx] = g_VBvolatile1(x,P,u,in)
% computes the probability of the learner choosing 'a=1' through
% probability matching
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
        X = [x(2);exp(x(3));x(7);exp(x(8))]; % [mu1;s1;mu2;s2]

gx = g_probability_matching( X,P,[],in );
% P(1): temperature/utility scaling
% P(2): preference bias towards choosing action a=1