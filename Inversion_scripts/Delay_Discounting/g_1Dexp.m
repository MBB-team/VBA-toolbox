function [gx] = g_1Dexp(x,P,u_t,in)


K = P(1);
r = exp(P(2));
beta = exp(P(3));
T = in.T;
V = in.V;

dU = ( (1-exp(-r*V(1,:)))./(1+K*T(1,:)) - (1-exp(-r*V(2,:)))./(1+K*T(2,:)) )./r;

gx = sig( beta*dU(:) );
end

function y=sig(x)
y = 1./(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;
end