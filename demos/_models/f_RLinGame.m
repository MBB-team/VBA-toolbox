function fx = f_RLinGame(x,P,u,in)
% Reinforcement-learner in a dyadic game
% function fx = f_RLinGame(x,P,u,in)
% An RL agent learns by trial and error. The RL agent updates its action
% values as follows:
% V(chosen action) = V(chosen action) + alpha*(feedback-V(chosen action))
% V(unchosen action) = V(unchosen action)
% In a dyadic game, the feedback is a function of both payers' actions. It
% is computed from the game's payoff table.
% NB: The payoff table is such that game(:,:,i) is the Ui(a,b), the payoff
% player i receives given that player 1 has played a and player 2 has
% played b.
% IN: 
%	- x: action values
%	- P: (invsigmoid-) learning rate
%	- u: u(1) = opponent's previous move, u(2) = agent's previous move
%   - in: structure:
%       u.game = 2x2x2 payoff table
%       u.player= agent's index (player's role)
% OUT:
%   - fx: updated action values
%   - dfdx/dfdP: gradients for VBA inversion

if VBA_isWeird (u) % e.g., 1st trial
    fx = x;
    return
end

% derive feedback from payoff table
player = in.player; % 1 or 2: role of the player
game = in.game; % payoff table
if player==1
    Payoff = game(:,:,1);
    P1 = 2-u(2); % player 1 = agent
    P2 = 2-u(1); % player 2 = opponent
elseif player==2
    Payoff = game(:,:,2);
    P1 = 2-u(1); % player 1 = opponent
    P2 = 2-u(2); % player 2 = agent
end
feedback = Payoff(P1,P2);

% update chosen action value given feedback
[fx,dfdx,dfdP] = f_Qlearn(x,P,[u(2);feedback],[]);
