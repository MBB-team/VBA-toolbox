function [gx] = g_Hampton(x,P,u,in)
% Hampton's "influence learning" model: observation function
% function [gx] = f_Hampton(x,P,u,in)
% Let P(o=1) be the agent's prediction of the other's next move, ie the
% probability that the other will pick the first alternative option. The
% "influence learner" bases her decision (a=1 or a=0) upon P(o=1), given
% the game payoff table.
% IN:
%   - x: hidden states:
%       x(1)= log-odds of P(o=1)
%   - P: parameters:
%       P(1)= (log-) temperature
%       P(2)= bias towards the first alternative option
%   - u: [useless]
%   - in: structure:
%       u.game = 2x2x2 payoff table
%       u.player= agent's index (player's role)
% OUT:
% - gx: proba that the agent will pick the first option, i.e. gx=P(a=1).

game = in.game; % game's payoff table
player = in.player; % agent's role
Po = VBA_sigmoid(x(1)); % P(o=1)
DV = fplayer(Po,exp(P(1)),player,game); % incentive for a=1
gx = VBA_sigmoid(DV+P(2)); % P(a=1) with bias