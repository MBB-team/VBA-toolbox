function [ gx ] = g_BSLinGame(x,P,u,in)
% observation function for BSL learner engaging in dyadic games
% [ gx ] = g_BSLinGame(x,P,u,in)
% In a game, a BSL learner bases her decision (a=1 or a=0) upon her
% prediction of her opponent's next move, given the game payoff table. NB:
% BSL tracks the probability that her opponent picks the first alternative
% option given the sequence of previous moves, ie P(ot=1|o{t-1K:t-1}).
% IN:
%   - x: sufficient statistics of log-odds of P(o=1):
%       x(1:2^K)= E[log-odds]
%       x((2^K)+1:2^(K+1))= log V[log-odds]
%   - P: observation param:
%       P(1) = log-temperature
%       P(2) = bias [optional]
%   - u: sequence of past actions:
%       u(1)= opponent's last move
%       u(2)= learner's last move
%       u(3:K+2) = sequence of K past opponent's moves
%   - in: structure:
%       u.game = 2x2x2 payoff table
%       u.player= agent's index (player's role)
% OUT:
%   - gx: proba that the agent will pick the first option, i.e. gx=P(a=1).

% Get the agent's prediction about her opponent's next move, ie P(o=1).
u(2) = []; % remove agent's last move
[Po] = g_BSL(x,[0;0],u,in);

% Make decision based upon the likely opponent's next move
DV = fplayer(Po,exp(P(1)),in.player,in.game); % incentive for a=1
gx = VBA_sigmoid(DV+P(2)); % P(a=1) with bias

