function [v] = v_discounting(x,P,u,in)
% delay discounting 2AFC observation function

t = u(1);
R = u(2);

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
        [tmp,i1] = min((in.gx-t).^2);
        [tmp,i2] = min((in.gy-R).^2);
        v = vec(in.bf(i1,i2,:))'*P;
end

