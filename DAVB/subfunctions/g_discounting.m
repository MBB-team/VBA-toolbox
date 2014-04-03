function [gx] = g_discounting(x,P,u,in)
% delay discounting 2AFC observation function

t = u(in.ind.t);
R = u(in.ind.R);

switch in.model
    case 'hyperbolic'
        k = exp(P(in.ind.logk));
        v = R./(1+k.*t);
    case 'exponential'
        k = exp(P(in.ind.logk));
        v = R.*exp(-k.*t);
    case 'linear'
        k = P(in.ind.logk);
        v = R - k*t;
end
dv = v(1) - v(2);
b = exp(-P(in.ind.logb));
gx = sigm(b.*dv);
