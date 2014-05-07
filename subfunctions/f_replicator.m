function [fx] = f_replicator(x,P,u,in)

% 1-evaluate fitness
h = feval(in.f_fitness,x,P,u,in);

% 2-evolve population
mh = mean(h);
f = h-mh;
fx = x + in.dt.*f.*x;



