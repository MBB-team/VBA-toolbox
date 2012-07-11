function [f,dfdx] = fitness_rsp(x,P,u,in)
% obtain fitness of hawk-dove game

c = P(1);

A = [   -c  1   -1
        -1  -c  1
        1   -1  -c   ];

% The malthusian fitness of player i is its payoff, averaged across the
% population. The chance of encountering any member of type j is
% proportional to its frequency wihin the population. This yields:
f = A*x;

dfdx = A;


