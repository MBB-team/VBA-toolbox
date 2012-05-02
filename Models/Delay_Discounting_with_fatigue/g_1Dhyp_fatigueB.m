function [gx] = g_1Dhyp_fatigueB(x,P,u_t,in)


K = exp(P(1));
beta1 = exp(P(2));
betaend = exp(P(3));

N = in.N; % total number of trials
beta = (betaend-beta1)/(N-1).*[1:N] + (beta1*N-betaend)/(N-1);

T = in.T;
V = in.V;

dU = V(1,:)./(1+K*T(1,:))-V(2,:)./(1+K*T(2,:));

gx = sigm( beta.*dU );

end

function y=sig(x)
y = 1./(1+exp(-x));
y(y<eps) = eps;
y(y>1-eps) = 1-eps;
end