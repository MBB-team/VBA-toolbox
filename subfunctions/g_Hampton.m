function [gx] = g_Hampton(x,P,u,in)
% Hampton's "influence learning" model: observation function
% function [gx] = f_Hampton(x,P,u,in)
% Hereafter, the learner is referred to as the "agent", whereas her
% opponent is referred to as the "other".
% IN:
%   - x: hidden states:
%       x(1)= estimated proba of other (needs to be sigmoid-mapped)
%       x(2)= FIXED temperature parameter (used for influence and choice)
%   - P: parameters:
%       P(1)= (log-) temperature
%       P(2)= bias towards the first alternative option
%   - u: [useless]
%   - in: cell-array:
%       u{1}= 2x2x2 payoff table
%       u{2}= agent's index (player's role)

game=in{1};
player=in{2};
Pi=sigmoid(x(1));
gx=sigmoid(fplayer(Pi,exp(P(1)),player,game)+P(2) );