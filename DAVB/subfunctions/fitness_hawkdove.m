function [f,dfdx] = fitness_hawkdove(x,P,u,in)
% obtain fitness of hawk-dove game

v = P(1);
c = P(2);

A = [   0.5*(v-c)   v
        0           0.5*v];

% The malthusian fitness of player i is its payoff, averaged across the
% population. The chance of encountering any member of type j is
% proportional to its frequency wihin the population. This yields:
f = A*x;

dfdx = A;


