function [gx] = g_kToM(x,P,u,in)
% wrapper around recursive k-ToM observation function
% function [gx] = g_kToM(x,P,u,in)
% IN:
%   - x: hidden states (see indexing in inG.indlev)
%   - P: observation param:
%       P(1) = (log-) temperature
%       P(2) = bias [optional]
%   - u: u(1) = opponent's previous move, u(2) = agent's previous move
%   - inG: input structure (see prepare_kToM.m)
% OUT:
%   - gx: proba that the agent will pick the first option, i.e. gx=P(y=1).
% [see ObsRecGen.m]

gx = ObsRecGen(x,P,u,in);