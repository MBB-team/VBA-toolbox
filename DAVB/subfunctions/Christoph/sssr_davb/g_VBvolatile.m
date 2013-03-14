function [gx] = g_VBvolatile(x,P,u,in)
% computes the probability of the learner choosing '1'.
% [gx] = g_VBvolatile(x,P,u,in)
% This function derives the action emission law from the expected next
% input to the learner. It is valid in an empirical context akin to a
% multiarmed bandit problem, whereby the learner is asked to bet on the
% outcome at each trial. This uses the states keeping trace of the previous
% posterior belief, at each trial, whose sufficient statistics are stored
% in the dummy states x(6) and x(7).
% Note: the action emission law accounts for uncertainty in the hidden
% states x2, when passed through the sigmoid mapping.
% IN:
%   - x: the previous posterior sufficient statistics. Here, only the dummy
%   states x(6) and x(7) serve in emitting a response at each trial.
%   - P: the response model parameters vector.
%   - u: the current input to the learner [not used here];
%   - in: [not used here].
% OUT:
%   - gx: the expected next input to the learner, according to its previous
%   updated belief about the hidden states.

ze = exp(P(1));
mu1hat = sgm(x(6),1);

gx = mu1hat.^ze./(mu1hat.^ze +(1-mu1hat).^ze);
