function  [fx] = f_wslsinGame(x,P,u,in)
% "win/stay - lose/switch" strategy in a dyadic game
% function  [fx,dfdx,dfdP] = f_wslsinGame(x,P,u,in)
% The "win/stay - lose/switch" strategy is encoded in terms of the
% evolution of pseudo q-values, which swap sign depending upon the
% feedback.
% In a dyadic game, the feedback is a function of both payers' actions. It
% is computed from the game's payoff table.
% NB: The payoff table is such that game(:,:,i) is the Ui(a,b), the payoff
% player i receives given that player 1 has played a and player 2 has
% played b.
% IN:
%   - x : pseudo q-values (1: stay, -1:switch)
%   - P : [useless]
%	- u: u(1) = opponent's previous move, u(2) = agent's previous move
%   - in: structure:
%       u.game = 2x2x2 payoff table
%       u.player= agent's index (player's role)
% OUT:
%   - fx: evolved pseudo q-values (2x1)

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

% apply "win/stay - lose/switch" heuristic given feedback
[fx] = f_wsls([],[],[u(2);feedback],[]);
