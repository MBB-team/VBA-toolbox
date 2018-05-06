function [ gx ] = ObsRecGen(x,P,u,in)
% observation function for k-ToM's bet about her opponent's next move
% [ gx ] = ObsRecGen(x,P,u,in)
% Marie Devaine wrote this function in November 2015 (comments: JD).
% A k-ToM learner bases her decision (a=1 or a=0) upon her prediction of
% her opponent's next move, given the game payoff table. When k>1, k-ToM
% maintains more than one such prediction, which depends upon the possible
% level of her opponent. Let P(o=1) be the probability that k-ToM's
% opponent will pick the first alternative option. Then:
% P(o=1) = sum_k P(o=1|k)*P(k)
% where P(o=1|k) is the probability that k-ToM's opponent will pick
% the first alternative option if he was a k-ToM, and P(k) is the
% probability that k-ToM's opponent is a k-ToM.
% IN:
%   - x: hidden states (see indexing in inG.indlev)
%   - P: observation param:
%       P(1) = (log-) temperature
%       P(2) = bias [optional]
%   - u: [useless]
%   - inG: input structure (see prepare_kToM.m)
% OUT:
%   - gx: proba that the agent will pick the first option, i.e. gx=P(y=1).

player = in.player; % 1 or 2: role of the player
ntotPar = in.npara; % only for k-ToM with k>0
level = in.lev; % depth of k-ToM's recursive beliefs
game = in.game; % payoff table
a = 0.36; % for E[s(x)] when x~n(mu,Sig)
indlev = in.indlev; % hidden-states indexing [see defIndlev.m]

% Get the agent's prediction about her opponent's next move, ie P(o=1).
if level==0 % 0-ToM
    
    mx = x(1); % E[log-odds of P(o=1)]
    Vx = exp(x(2)); % V[log-odds of P(o=1)]
    Po = VBA_sigmoid(mx/(sqrt(1+a*Vx))); % P(o=1)
    
else
    
    % Get P(k'). Note: if the agent is k-ToM, then, by definition, she
    % considers that her opponent's sophistication is k' < k. In addition,
    % there is a constraint of normalization, ie sum_k' P(k') = 1. Thus,
    % one only needs to keep track of k'-1 probabilities (the last one is,
    % by construction, 1-sum_k' P(k')).
    Pk = VBA_sigmoid(x(1:(level-1))); % P(k'), with k'=0,...,k-1
    Pk = [Pk;max(0,1-sum(Pk))]; % insert last P(k'=k-1)
    
    % Get P(o=1|k'). Note: the agent's prediction P(o=1|k') depends upon
    % her estimate of her opponent's parameters (learning rate, tmperature,
    % bias...). Uncertainty Re: these parameters eventually results in
    % blurring her prediction.
    % Note: hidden states encode x(theta), the log-odds of P(o=1|k',theta)
    % evaluated at the agent's estimate of theta (mu). In addition, they
    % encode the gradient of x wrt to theta (dx/dtheta), and V[theta]. This
    % then serves to derive P(o=1|k') as follows:
    % P(o=1|k') = E[sigm(x(theta))]
    %           = sigm(E[x(theta)]/sqrt(1+a*V[x(theta)])
    %           = sigm(E[x(theta)]/sqrt(1+a*V[theta]*(dx/dtheta)^2)   
    f = zeros(level,1); % E[x(theta)]
    Vx = zeros(level,1); % V[x(theta)]
    for j=1:level % loop over possible opponent's levels  (k'=j-1)
        f(j) = x(indlev(j).f); % E[x(theta)|k'=j-1]
        df = x(indlev(j).df); % d[x(theta)]/dtheta for k'=j-1
        Sig = exp(x(indlev(j).Par(2:2:2*ntotPar))); % V[theta|k'=j-1]
        Vx(j) = sum(Sig.*df.^2); % V[x(theta)|k'=j-1]
    end
    Es = VBA_sigmoid(f./sqrt(1+a*Vx)); % E[sigm(x(theta))]
    
    % Get P(o=1) = sum_k P(o=1|k')*P(k')
    Po = Pk'*Es; % k-ToM's belief about her opponent's next move
    
end

% Make decision based upon P(o=1)
DV = fplayer(Po,exp(P(1)),player,game); % incentive for a=1
if length(P)==1
    gx = VBA_sigmoid(DV); % P(a=1)
else % P(2) = bias
    gx = VBA_sigmoid(DV+P(2)); % P(a=1) with bias
end

