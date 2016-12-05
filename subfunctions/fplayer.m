function DV=fplayer(P,beta,player,game)
% derives the decision variable for a player relying on a 2x2 game payoff
% function DV=fplayer(P,beta,player,game)
% IN:
%   - P: probability that the opponent's action will be the first option
%   - beta: behavioural temperature
%   - player: player index (1 or 2)
%   - game: 2x2x2 payoff matrix
% OUT:
%   - DV: decision variable

% if player==2
%     Payoff= game(:,:,2);
%      DV=((Payoff(1,1)+Payoff(2,2))*P+(1-P)*Payoff(2,1)-Payoff(1,2)*P-Payoff(2,2))/beta;  
% elseif player==1
%     Payoff= game(:,:,1);
%     DV=((Payoff(1,1)+Payoff(2,2))*P+(1-P)*Payoff(1,2)-Payoff(2,1)*P-Payoff(2,2))/beta;  
% end

if player==2
    Payoff= game(:,:,2);
    DV = P*(Payoff(1,1)-Payoff(1,2))+(1-P)*(Payoff(2,1)-Payoff(2,2));
elseif player==1
    Payoff= game(:,:,1);
    DV = P*(Payoff(1,1)-Payoff(2,1))+(1-P)*(Payoff(1,2)-Payoff(2,2));
end
DV = DV/beta;
