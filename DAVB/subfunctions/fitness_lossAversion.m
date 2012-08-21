function f = fitness_lossAversion(x,P,u,in)
% obtain fitness of loss-averse types playing stag-hunt

% Loss-aversion trait is defined on a 1D continuous interval, which is the
% range of a parameter controlling the gain/loss asymmetry of the utility
% function of agents (see function lossAverseUtile.m). When this parameter
% is equal to zero (resp., one), agents are insensitive to gains (resp.,
% losses).
% Let lambda be the parameter controlling the gain/loss asymmetry, and
% x(lambda) be the frequency of this trait within the population. This demo
% uses standard replicator dynamics of evolutionary game theory to derive
% the adaptive fitness of the trait, ie:
%
% dx_i/dt = x_i ( f_i(x) + sum_j{x_j*f_j(x)} )
%
% where f_i(x) is the malthusian fitness of trait i (ie of a given lambda).
% This fitness is obtained according to the expected outcome of a game
% (here: "stag-hunt"), for each type of agent or trait.
% In our case, lambda changes the probability to cooperate, through a
% probabilistic softmax action emission law.
%
% The outcome table of the "stag-hunt" game is as follows:
% A = [ X(1)-c  X(3)-c
%       X(2)-c  X(2)-c  ]
% where the first line of A applies when the player cooperates, and the
% second line applies when she defects (columns pertain to the opponent's
% behaviour). Note that the payoff of defecting is independent of the
% opponent's behaviour.
% There is an effort cost (c) attached to any action. This cost does not
% affect the action emission law, since the difference in values is
% invariant to c. However, it turns out to be critical for deriving the
% evolutionary fitness of loss aversion.


n = size(x,1);

% form outcome table for player 1:
A = [   in.X(1) in.X(3)
        in.X(2) in.X(2) ];
A = A - in.mc + in.sc.*randn(1,1);

% get proba of choosing stag for each population
p = zeros(n,1);
for i=1:n
    lambda = in.lambda(i);
    pA = lossAverseUtile(A,lambda)*ones(2,1)/2;
    du = pA(1)-pA(2);
    p(i) = 1./(1+exp(-du));
end

% get mean outcome for each pair of population
O = zeros(n,n);
for i=1:n
    % get proba, for trait i, of playing cooperate and defect (p,1-p)
    Pi = [p(i);1-p(i)];
    for j=1:n
        % get proba, for trait j, of playing cooperate and defect (p,1-p)
        Pj = [p(j);1-p(j)];
        % get expecte outcome of i, when playing against against j
        O(i,j) = Pi'*A*Pj;
    end
end
% The malthusian fitness of player i is its payoff, averaged across the
% population. The chance of encountering any member of type j is
% proportional to its frequency wihin the population. This yields:
f = O*x;



