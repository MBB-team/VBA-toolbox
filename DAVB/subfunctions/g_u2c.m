function gx = g_u2c(x,P,u,in)
% generates the probability of picking 1 item among 2 from parameterized
% utility function

dv = P(in.ind(:,1)) - P(in.ind(:,2)); % relative value of chosen item
b = exp(P(in.temp)); % behavioural temperature
gx = sig(dv/b); % probability of picking the first item

function s= sig(x)
s = 1./(1+exp(-x));
s(s<1e-3) = 1e-3;
s(s>1-1e-3) = 1-1e-3;