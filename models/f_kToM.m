function [fx] = f_kToM(x,P,u,inF)
% wrapper around recursive k-ToM evolution function
% function [fx] = f_kToM(x,P,u,inF)
% A k-ToM agent with k>0 learns the unknown parameters (theta) of her
% opponent's learning and decision rules. Typically, theta includes her
% opponent's prior volatility about herself (which controls the learning
% rate), the behavioural temperature and the bias. In addition, when k>1,
% k-ToM also learns her opponent's sophistication level (k').
% Note that k-ToM assumes that her opponent's parameters can drift over
% trials. How much they can change is controlled by her prior volatility
% about her opponent's parameters. In addition, the agent may partially
% "forget" about her opponent's sophistication level from a trial to the
% next. This effectively dilutes the belief P(k') towards the corresponding
% max-entropic distribution (only if inF.diluteP=1).
% IN:
%   - x: hidden states (see indexing in inF.indlev)
%   - P: evolution params:
%   P(1)= agent's prior opponent's (log-) volatility on opponent's params
%   P(2)= agent's (invsigmoid-) dilution coefficient
%   - u: u(1) = opponent's previous move, u(2) = agent's previous move
%   - inF: input structure (see prepare_kToM.m)
% OUT:
%   - fx: updated hidden states
% [see RecToMfunction.m]

if VBA_isWeird (u) % e.g., 1st trial
    fx = x;
    return
end

[fx,indlev] = RecToMfunction(x,P,u,inF);
% NB: output indlev is used recursively by RecToMfunction.m