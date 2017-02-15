function [gx] = g_kToM(x,P,u,in)
% wrapper around recursive k-ToM observation function
% function [gx] = g_kToM(x,P,u,in)
% A k-ToM learner bases her decision (a=1 or a=0) upon her prediction of
% her opponent's next move, given the game payoff table.
% IN:
%   - x: hidden states (see indexing in inG.indlev)
%   - P: observation param:
%       P(1) = (log-) temperature
%       P(2) = bias [optional]
%   - u: [useless]
%   - inG: input structure (see prepare_kToM.m)
% OUT:
%   - gx: proba that the agent will pick the first option, i.e. gx=P(a=1).
% [see ObsRecGen.m]

gx = ObsRecGen(x,P,u,in); % P(a=1)