function gx = g_u2c(x,P,u,in)
% generates the probability of picking 1 item among 2 from parameterized
% utility function

% alternatives
ac = u(in.ic); % item 1
au = u(in.iu); % item 2
if ~isempty(in.effort)
    k = exp(P(in.effort.P)); % effort cost per force unit
    c1 = -k*u(in.effort.ic); % efort cost of item 1
    c2 = -k*u(in.effort.iu); % efort cost of item 2
else
    c1 = 0;
    c2 = 0;
end
dv = P(ac) + c1 - P(au) - c2; % relative value of item 1
b = exp(P(in.temp)); % behavioural temperature
gx = sig(dv/b); % probability of picking the first item

function s= sig(x)
s = 1./(1+exp(-x));
s(s<1e-3) = 1e-3;
s(s>1-1e-3) = 1-1e-3;