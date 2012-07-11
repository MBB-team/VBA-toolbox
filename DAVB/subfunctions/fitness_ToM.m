function f = fitness_ToM(x,P,u,in)
% obtain fitness of loss-averse types playing stag-hunt

lambda = P(1);
tau = P(2);

[tmp,it] = min((in.tau-tau).^2);
A = lambda*in.Q(:,:,it,in.coop_game) + (1-lambda)*in.Q(:,:,it,2);

% The malthusian fitness of player i is its payoff, averaged across the
% population. The chance of encountering any member of type j is
% proportional to its frequency wihin the population. This yields:
f = A*x;





