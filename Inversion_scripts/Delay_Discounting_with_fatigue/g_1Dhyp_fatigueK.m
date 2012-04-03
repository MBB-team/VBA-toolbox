function [gx] = g_1Dhyp_fatigueK(x,P,u_t,in)


K1 = exp(P(1));
Kend = exp(P(2));
beta = exp(P(3));

N = in.N; % total number of trials
K = (Kend-K1)/(N-1).*[1:N] + (K1*N-Kend)/(N-1);

T = in.T;
V = in.V;

dU = V(1,:)./(1+K.*T(1,:))-V(2,:)./(1+K.*T(2,:));

gx = sig( beta*dU(:) );
end

function y=sig(x)
y = 1./(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;
end