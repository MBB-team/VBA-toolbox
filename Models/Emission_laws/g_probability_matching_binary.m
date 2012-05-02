function [gx] = g_probability_matching(x,P,u,in)
% computes the probability of the learner choosing '1'.
% "A bayesian foundation for individual learning under uncertainty"
% IN:
%   - x: the previous posterior sufficient statistics. Here, only the dummy
%   states x(6) and x(7) serve in emitting a response at each trial.
%   - P: the response model parameters vector.
%   - u: the current input to the learner [not used here];
%   - in: [not used here].
% OUT:
%   - gx: the expected next input to the learner, according to its previous
%   updated belief about the hidden states.
beta = exp(P(1));
gx = x(6)^beta/(x(6)^beta+(1-x(6))^beta);


