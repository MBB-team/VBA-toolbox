function [ gx ] = g_HGFinGame(x,P,u,in)
% observation function for HGF learner engaging in dyadic games
% [ gx ] = g_HGFinGame(x,P,u,in)
% In a game, a HGF learner bases her decision (a=1 or a=0) upon her
% prediction of her opponent's next move, given the game payoff table. NB:
% HGF tracks the probability that her opponent picks the first alternative
% option, ie P(o=1).
% IN:
%   - x: posterior sufficient statistics:
%   x(1)= last o
%   x(2)= E[log-odds of P(o=1)]
%   x(3)= log V[log-odds of P(o=1)]
%   x(4)= E[log-volatility]
%   x(5)= log V[log-volatility]
%   - P: observation param:
%       P(1) = (log-) temperature
%       P(2) = bias [optional]
%   - u: [useless]
%   - in: structure:
%       u.game = 2x2x2 payoff table
%       u.player= agent's index (player's role)
% OUT:
%   - gx: proba that the agent will pick the first option, i.e. gx=P(a=1).

player = in.player; % 1 or 2: role of the player
game = in.game; % payoff table
a = 0.36; % for E[s(x)] when x~n(mu,Sig)

% Get the agent's prediction about her opponent's next move, ie P(o=1).
mx = x(2); % E[log-odds of P(o=1)]
Vx = exp(x(3)); % V[log-odds of P(o=1)]
Po = VBA_sigmoid(mx/(sqrt(1+a*Vx))); % P(o=1)

% Make decision based upon the likely opponent's next move
DV = fplayer(Po,exp(P(1)),player,game); % incentive for a=1
gx = VBA_sigmoid(DV+P(2)); % P(a=1) with bias


