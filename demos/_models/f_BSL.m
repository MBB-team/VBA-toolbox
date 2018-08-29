function [fx] = f_BSL(x,theta,u,in)
% evolution function for a Bayesian sequence-learner (BSL)
% function [fx] = f_BSL(x,theta,u,in)
% BSL is simply tracking the log-odds of P(y_t=1|y_{t-1}), where y is a
% binary outcome. This variable is updated according to a Laplace-Kalman
% filter, yielding 2 sufficient statistics (m and V) per combination of
% past outcome. BSL can learn sequences of arbitrary depth (K). For
% example, if K=1, then BSL tracks 2 probabilities, namely:
% P(y_t=1|y_{t-1}=1) and P(y_t=1|y_{t-1}=0). More generally, BSL tracks 2^K
% probabilities. In this scheme, the only evolution param (theta) is BSL's
% prior volatity about the log-odds.
% Note: unsampled sequences will eventually be "forgotten", since the
% prediction step in the Laplace-Kalman update will dilute any previously
% sampled evidence.
% IN:
%   - x: sufficient statistics of log-odds of P(y=1):
%       x(1:2^K)= E[log-odds]
%       x((2^K)+1:2^(K+1))= log V[log-odds]
%   - theta: BSL's prior volatity
%   - u: u(1)= current outcome, u(2:K+1)= sequence of K past outcomes
%   - in: depth of sequence learning
% OUT:
%   - fx: updated sufficient statistics of log-odds of P(y=1)

if VBA_isWeird (u) % e.g., 1st trial
    fx = x;
    return
end

% -- deal with superceding competition with other generative models --
% [See, e.g., f_metaTom.m]
try
    w = inF.metaweight;
catch
    w = 1;
end

% -- learning rule --
K = in.K; % sequence depth
y = u(1); % current outcome
yb = u(2:K+1); % previous outcomes
m0 = x(1:2^K); % current E[log-odds]
V0 = exp(x((2^K)+1:2^(K+1))); % current V[log-odds]
volatility = exp(theta(1)); % positivity constraint
m = m0; % by default: no change in E[log-odds]
if K >0
    indSeq = bin2dec(num2str(yb'))+1; % index of sequence of previous outcomes
else
    indSeq = 1;
end
V = V0 + volatility; % by default: inflated V[log-odds]
p0 = VBA_sigmoid(m0(indSeq)); % current estimate of P(y=1)
V(indSeq) = 1./((1./V(indSeq))+ w*p0*(1-p0)); % updated V[log-odds]
m(indSeq) = m0(indSeq) + w*V(indSeq)*(y-p0); % updated E[log-odds] (Laplace-Kalman update rule)
fx = [VBA_sigmoid(VBA_sigmoid(m),'inverse',true);log(V)]; % wrap-up


 