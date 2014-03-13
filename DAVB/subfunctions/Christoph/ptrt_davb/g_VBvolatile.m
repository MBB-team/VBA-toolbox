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
%   - u: the current input to the learner.
%   - in: further quantities handed to the function.
% OUT:
%   - gx: the expected next input to the learner, according to its previous
%   updated belief about the hidden states.

respmod = in.respmod;

ze1 = exp(P(1));
ze2 = exp(P(2));

mu2hat = x(6);
signmu2hat = sign(mu2hat);
mu1hat = sgm(mu2hat,1);
pi1hat = 1/(mu1hat*(1-mu1hat));

switch respmod
  case 'precision'
    gx = ze1+ze2./(1+exp(-signmu2hat.*(2.*u-1).*(pi1hat-4)));
  case 'surprise'
    gx = ze1+ze2./(1-u.*log(mu1hat)-(1-u).*log(1-mu1hat));
  case 'belief'
    gx = ze1+ze2.*mu1hat.^u.*(1-mu1hat).^(1-u);
  otherwise
    error(['Invalid or missing response model specification: ', respmod]);
end
