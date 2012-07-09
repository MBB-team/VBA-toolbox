function [gx] = g_1Dhyp(x_t,P,u_t,in)
%%% Binary choice, hyperbolic delay discounting + softmax decision
% INPUT
% - x_t : []
% - P : 2*1, discount factor, inverse temperature in softmax (both
% exponential mapping)
% - u_t : []
% - in :
%   - V : 2*Ntrials, 'objective' value of proposed alternatives
%   - T : 2*Ntrials, time of delivery of the proposed alternatives
% OUTPUT
% - gx : Scalar, P(a=a1|x_t), probability of choice 1

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