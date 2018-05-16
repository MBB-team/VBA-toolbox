function DV = fplayer(P,beta,player,game)
% derives the player's decision variable given a 2x2x2 payoff table
% function DV = fplayer(P,beta,player,game)
% The payoff table is such that game(:,:,i) is the Ui(a,b), the payoff
% player i receives given that player 1 has played a and player 2 has
% played b. The agent holds a belief about her opponent's next move, in
% terms of the probability p=P(o=1) that her opponent will play the first
% alternative action. Her decision variable DV is given by:
%  DV = p*(Ui(a=1,b=1)-Ui(a=0,b=1)) + (1-p)*(Ui(a=1,b=0)-Ui(a=0,b=0))
% where actions have been given binary values by convention.
% IN:
%   - P: probability that the opponent's action will be the first option
%   - beta: behavioural temperature
%   - player: player index (1 or 2)
%   - game: 2x2x2 payoff matrix
% OUT:
%   - DV: decision variable (incentive for chosing option 1)

if player==2
    Payoff = game(:,:,2);
    DV = P*(Payoff(1,1)-Payoff(1,2))+(1-P)*(Payoff(2,1)-Payoff(2,2));
elseif player==1
    Payoff = game(:,:,1);
    DV = P*(Payoff(1,1)-Payoff(2,1))+(1-P)*(Payoff(1,2)-Payoff(2,2));
end
DV = DV./beta;

