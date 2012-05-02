function [gx] = g_Mathys_multinomial_inv(x,P,u,in)
% computes the probability of the learner choosing '1'.
% "A bayesian foundation for individual learning under uncertainty"
% IN:
%   - x: the previous posterior sufficient statistics. Here, only the dummy
%   states x(6) and x(7) serve in emitting a response at each trial.
%   - P: the response model parameters vector.
%   - u: the current input to the learner [not used here];
%   - in: [not used here].
% OUT:
%   - gx: the probabilities of choosing the different actions
%   This is a vector of dimension Na (Na = number of alternatives)

v = g_Mathys_multinomial_sim(x,P,u,in);
gx = v(u(1));


