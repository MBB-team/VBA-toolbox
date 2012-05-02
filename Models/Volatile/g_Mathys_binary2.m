function [gx] = g_Mathys_binary2(x,P,u,in)
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

Na = size(x,1)/7;
b = exp(P(1));
v = x(6:7:end,1); % vector of proba of outcomes for each alternatives
gx = exp( b*v(u(1)) )/sum( b*exp(v) );

%sgm(x(6),1) 

