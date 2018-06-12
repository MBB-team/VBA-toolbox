function [gx] = g_VBvolatile0(x,P,u,in)
% computes the probability of the learner choosing 'a=1'.
% [gx] = g_VBvolatile0(x,P,u,in)
% This function derives the action emission law from the current belief of
% the learner. It is valid in an empirical context akin to a multiarmed
% bandit problem, whereby the learner is asked to bet on the outcome at
% each trial.
% Note: the state-space includes the first- and second- order moments of
% the log-odds probability that each action yield a reward, ie P(R|a).
% IN:
%   - x: the previous posterior sufficient statistics.
%   - P: the response model parameters vector:
%   P(1): temperature/utility scaling
%   P(2): preference bias towards choosing action a=1
%   - u: the current input to the learner.
%   - in: further quantities handed to the function.
% OUT:
%   - gx: the expected next input to the learner, according to its previous
%   updated belief about the hidden states.


switch in.respmod
    case 'taylor' % neglect uncertainty on P(R|a=1)
        x1 = x(2); % E[log-odds P(R|a=1)]
        x0 = x(7); % E[log-odds P(R|a=0)]
    case 'fixedForm' % account for uncertainty on P(R|a=1)
        a = 0.368;
        x1 = x(2)./sqrt(1+a*exp(x(3)));
        x0 = x(7)./sqrt(1+a*exp(x(8)));
    otherwise
        error(['Invalid or missing response model specification: ', respmod]);
end

p1 = VBA_sigmoid(x1); % E[P(R|a=1)]
p0 = VBA_sigmoid(x0); % E[P(R|a=0)]
gx = VBA_sigmoid((p1-p0)*exp(P(1))+P(2)); % P(a=1)

