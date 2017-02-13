function [fx] = f_kToM(x,P,u,inF)
% wrapper around recursive k-ToM evolution function
% function [fx] = f_kToM(x,P,u,inF)
% IN:
%   - x: hidden states (see indexing in inF.indlev)
%   - P: evolution param = prior opponent's volatility
%   - u: u(1) = opponent's previous move, u(2) = agent's previous move
%   - inF: input structure (see prepare_kToM.m)
% OUT:
%   - fx: updated hidden states
% [see RecToMfunction.m]

[fx,indlev] =RecToMfunction(x,P,u,inF);