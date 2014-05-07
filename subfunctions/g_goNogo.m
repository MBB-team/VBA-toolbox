function [gx] = g_goNogo(x,P,u,in)
% generates the probability of 'go' choice: P(go) = sig(V/temperature+bias)
% function [gx] = g_goNogo(x,P,u,in)
% IN:
%   - x: value of the 'go' option
%   - P: log-temperature and bias
%   - u: [useless]
%   - in: [useless]
% OUT:
%   - gx: P(go) = sig(V/temperature+bias)

gx = 1./(1+exp(P-x));