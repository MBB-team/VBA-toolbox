function gx = g_u2p(x,P,u,in)
% generates the probability of picking 1 item among 2 from parameterized
% utility function

mu = x(1:in.n); % expectation (E[x] = mu)
dv = mu(u(in.iu(1))) - mu(u(in.iu(2))); % relative value of item 1
b = exp(P(in.temp)); % behavioural temperature
gx = VBA_sigmoid(dv/b, 'finite', 1e-3); % probability of picking the first item
