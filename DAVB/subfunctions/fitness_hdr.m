function [f,dfdx] = fitness_hawkdove(x,P,u,in)
% obtain fitness of hawk-dove(retaliator game

c = P(1);

A = [   -1  2   -1
        0   1   1-c
        -1  1+c 1  ];
    
% The malthusian fitness of player i is its payoff, averaged across the
% population. The chance of encountering any member of type j is
% proportional to its frequency wihin the population. This yields:
f = A*x;

dfdx = A;


