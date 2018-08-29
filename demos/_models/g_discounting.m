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
    case 'basis'
        for i=1:2
            [tmp,i1(i)] = min((in.grid1-t(i)).^2);
            [tmp,i2(i)] = min((in.grid2-R(i)).^2);
            v(i) = VBA_vec(in.bf(i1(i),i2(i),:))'*P;
        end
end
dv = v(1) - v(2);
b = exp(-P(in.ind.logb));
gx = VBA_sigmoid(b.*dv);
