function [fx] = f_Hampton(x,P,u,in)
% Hampton's "influence learning" model: evolution function (Hide and Seek)
% function [fx] = f_Hampton(x,P,u,in)
% Hereafter, the learner is referred to as the "agent", whereas her
% opponent is referred to as the "other".
% Let p=P(o=1) be the agent's prediction of the other's next move, ie the
% probability that the other will pick the first alternative option. The
% "influence learning" rule can be written as follows:
% p = p0 + eta*(o-p0) - lambda*k1*p0*(1-p0)*(a-q0)
% where o is the other's last move, a is the agent's last move, eta is the
% weight of the agent's prediction error, lambda is the weight of other's
% prediction error, k1 is derived from the game's payoff table. Here, q0 is
% the opponent's belief about the agent last move, and it is derived from
% p0, the the game's payoff table and the opponent's temperature.
% IN:
%   - x: hidden states:
%       x(1)= log-odds of P(o=1)
%   - P: learning parameters:
%       P(1)= (invsigmoid-) weight of agent's prediction error (invsigmoid-)
%       P(2)= (invsigmoid-) weight of other's prediction error
%       P(3)= (log-) other's temperature
%   - u: u(1) = opponent's previous move, u(2) = agent's previous move
%   - in: structure:
%       u.game = 2x2x2 payoff table
%       u.player= agent's index (player's role)
% OUT:
%   - fx: updated hidden states

if VBA_isWeird (u) % e.g., 1st trial
    fx = x;
    return;
end

o = u(1); % opponent's last choice (o)
a = u(2); % agent's last choice (a)
p0 = VBA_sigmoid(x(1)); % previous estimate of P(o=1)
eta = VBA_sigmoid(P(1)); % weight of PE1
lambda = VBA_sigmoid(P(2)); % weight of PE2
beta = exp(P(3)); % opponent's temperature

% derive first-order prediction error
PE1 = o-p0;

% derive second-order prediction error
game = in.game; % game's payoff table
player = in.player; % agent's role (opponent's role = 3-player)
if player==2
    Payoff = game(:,:,1);
    k1 = Payoff(1,1)-Payoff(2,1)-Payoff(1,2)+Payoff(2,2);
    k2 = Payoff(2,1)-Payoff(2,2);
elseif player==1
    Payoff = game(:,:,2);
    k1 = Payoff(1,1)-Payoff(2,1)-Payoff(1,2)+Payoff(2,2);
    k2 = Payoff(1,2)-Payoff(2,2);
end
q0 = (beta.*x(1) - k2)./k1;
PE2 = a-q0;

% "influence" learning rule
p = p0 + eta*PE1 + lambda*k1*p0*(1-p0)*PE2; % P(o=1)
p = max([min([p,1]),0]); % bound p between 0 and 1
fx = VBA_sigmoid(p,'inverse',true); % for numerical reasons

