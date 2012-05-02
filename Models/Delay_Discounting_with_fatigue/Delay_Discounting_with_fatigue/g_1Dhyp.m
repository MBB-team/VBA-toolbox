function [gx] = g_1Dhyp(x,P,u_t,in)


K = exp(P(1));
beta = exp(P(2));

T = in.T;
V = in.V;

dU = V(1,:)./(1+K*T(1,:))-V(2,:)./(1+K*T(2,:));

gx = sig( beta*dU(:) );
end

function y=sig(x)
y = 1./(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;
end