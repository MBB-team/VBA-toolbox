function [fx] = log_replicator(xlo,P,u,in)
% evaluates the replicator mapping of EGT (in log-odds space)
% function [fx] = log_replicator(xlo,P,u,in)
% IN:
% - xlo: log-odds of competitive traits within the population
% - P: evolution parameters
% - u: input to the system
% - in: optional structure containing the handle of the function that
% evaluates the fitness of each trait (in.f_fitness)
% OUT:
% - fx: the evolved log-odds
%
% Let x_k be the frequency of level-k players within the population.
% Standard replicator dynamics of evolutionary game theory write:
%
% dx_k/dt = x_k ( f_k(x) + sum_j{x_j*f_j(x)} )
%
% where f_k(x) is the malthusian fitness of the k^th trait.
% This equation however, does not ensure proper normalization of trait
% frequencies, ie: x_k > 0 and sum_j{x_j} = 1.
% We therefore define x_k to be the softmax mapping of log-odds xlo, ie:
%
% x_k = s_k(x_lo)
%     = exp(xlo_k)./sum_j{exp(xlo_j)}
%
% This means we can derive the dynamics of frequency traits under proper
% normalization constraints by integrating the evolution function of
% log-odds in time, and applying the above softmax mapping:
%
% dxlo_k/dt = [ds_k/dx]^-1*dx_k/dt = f_k(s(xlo))
% x(t)      = s(x_lo(t))
% 
% This function furnishes the Euler approximate step for the above
% evolution function in log-odds space, given the malthusian fitness of
% each trait witin the population.

% 0-sigmoid-transform log-odds
[x] = g_odds(xlo,P,u,in);

% 1-evaluate fitness
h = feval(in.f_fitness,x,P,u,in);

% % 2-evolve population
% mh = mean(h);
% f = h-mh;
% 
% % 3- apply inverse Jacobian of softmax mapping
% n = size(x,1);
% I = eye(n);
% J = I- ones(n,1)*x';
% iJ = pinv(J);
% f = iJ*f;

f = h;

% 4- apply Euler discretization step
fx = xlo + in.dt.*f;

% 5- remove the mean for numerical stability
fx = fx - mean(fx);



