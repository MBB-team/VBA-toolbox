function gx = g_u2c(x,P,u,in)
% generates the probability of picking 1 item among 2 from parameterized
% utility function

dv = P(in.ind(:,1)) - P(in.ind(:,2)); % relative value of chosen item
b = exp(P(in.temp)); % behavioural temperature
gx = VBA_sigmoid(dv/b, 'finite', 1e-3); % probability of picking the first item