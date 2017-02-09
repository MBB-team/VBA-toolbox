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
%       x(1)= log-odds of p
%   - P: parameters:
%       P(1)= (log-) weight of agent's prediction error (~=learning rate)
%       P(2)= (log-) weight of other's prediction error
%       P(3)= other's (log-) temperature
%   - u: inputs:
%       u(1)= y2 (other's last choice)
%       u(2)= y1 (agent's last choice)
%   - in: cell-array:
%       u{1}= 2x2x2 payoff table
%       u{2}= agent's index (player's role)

y2=u(1);% choice Other pmayer
y1=u(2); % choice Player 1
p0=sigmoid(x(1));
eta = exp(P(1));
lambda = exp(P(2));
beta = exp(P(3));
PE1 = y2-p0;

game = in{1};
player = in{2};
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
PE2 = y1-q0;

p = p0 + eta*PE1 + lambda*k1*p0*(1-p0)*PE2;

% % q_2=max(min( .5*(1-x(1)*exp(P(3))),1),0);
% q_2 = .5*(1-x(1)*exp(P(3)));
% p=p0+ exp(P(1))*(y2-p0)-2*exp(P(2))*p0*(1-p0)*(y1-q_2)
% pause

fx = min([max([p,0]),1]); % bound p between 0 and 1

