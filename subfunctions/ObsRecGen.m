function [ gx ] = ObsRecGen(x,P,u,in)
% observation function for k-ToM's bet about her opponent's next move
% [ gx ] = ObsRecGen(x,P,u,in)
% Marie Devaine wrote this function in November 2015 (comments: JD).
% A k-ToM learner bases her decision upon her prediction of her opponent's
% next move, given the game payoff table. When k>1, k-ToM maintains more
% than one such prediction, which depends upon the possible level of her
% opponent. Let P(o=1) be the probability that k-ToM's opponent will pick
% the first alternative option. Then:
% P(o=1) = sum_k P(o=1|k)*P(k)
% where P(o=1|k) is the probability that k-ToM's opponent will pick
% the first alternative option if he was a k-ToM, and P(k) is the
% probability that k-ToM's opponent is a k-ToM.
% IN:
%   - x: hidden states (see indexing in inG.indlev)
%   - P: observation param:
%       P(1) = (log-) temperature
%       P(2) = bias [optional]
%   - u: u(1) = opponent's previous move, u(2) = agent's previous move
%   - inG: input structure (see prepare_kToM.m)
% OUT:
%   - gx: proba that the agent will pick the first option, i.e. gx=P(y=1).

player = in.player; % 1 or 2: role of the player
ntotPara = in.npara; % only for k-ToM with k>0
lev = in.lev; % depth of k-ToM's recursive beliefs
game = in.game; % payoff table
a = 0.36; % for E[s(x)] when x~n(mu,Sig)
indlev = in.indlev; % hidden-states indexing [see defIndlev.m]

% Get the agent's prediction about her opponent's next move, ie P(o=1).
if lev==0 % 0-ToM
    
    muT = x(1); % E[log-odds of P(o=1)]
    SigT = exp(x(2)); % V[log-odds of P(o=1)]
    Pi = sigmoid( muT/(sqrt(1+a*SigT))); % P(o=1)
    
else
    
    % Get P(k). Note: if the agent is k'-ToM, then, by definition, she
    % considers that her opponent's sophistication is k <= k'-1. In
    % addition, there is a constraint of normalization, ie sum_k P(k) = 1.
    % Thus, one only needs to keep track of k'-1 probabilities (the last
    % one is, by construction, 1-sum_k P(k)).
    vecp = vec(sigmoid(x(1:(lev-1)))); % P(k), with k=0,...,k'-2
    vecp = [vecp;max(0,1-sum(vecp))]; % insert last P(k=k'-1)
    
    % Get P(o=1|k). Note: the agent's prediction P(o=1|k) depends upon her
    % estimate of her opponent's parameters (learning rate, tmperature,
    % bias...). Uncertainty Re: these parameters eventually results in
    % blurring her prediction, because:
    % P(o=1|k) = Integral_theta { P(o=1|k,theta) dtheta }
    % where theta are her opponent's parameters.
    % Here, hidden states store x(theta) = log-odds of P(o=1|k,theta),
    % its derivative wrt theta, and V[theta]. 
    Vx = zeros(lev,1); % E[x(theta)]
    vecPi = zeros(lev,1);
    for j=1:lev % loop over the opponent's possible levels (k=j-1)
        f = x(indlev(j).f); % E[x(theta)|k=j-1]
        df = x(indlev(j).df)'; % d[x(theta)]/dtheta for k=j-1
        for i=1:ntotPara % loop over the opponent's params
            Sig = exp(x(indlev(j).Par(2*i))); % V[theta|k=j-1]
            Vx(j) = Vx(j)+Sig*df(i)^2; % V[x(theta)|k=j-1]
        end
        vecPi(j)=sigmoid(f/sqrt(1+a*Vx(j))); % P(o=1|k) = E[s(f)]
    end
    
    % Get P(o=1) = sum_k P(o=1|k)*P(k)
    Pi = vecp'*vecPi;
    
end

% Make decision based upon P(o=1)
DV = fplayer(Pi,exp(P(1)),player,game); % incentive for a=1
if length(P)==1
    gx = sigmoid(DV); % P(a=1)
else % P(2) = bias
    gx = sigmoid(DV+P(2)); % P(a=1) with bias
end

