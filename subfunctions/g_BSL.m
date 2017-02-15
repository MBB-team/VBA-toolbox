function [gx] = g_BSL(x,phi,u,in)
% observation function for a Bayesian sequence-learner (BSL)
% function [gx] = g_BSL(x,theta,u,in)
% BSL is guessing the next outcome y_t based upon P(y_t=1|y_{t-1}), ie her
% bet actually depends upon the past sequence of outcomes.
% IN:
%   - x: sufficient statistics of log-odds of P(o=1):
%       x(1:2^K)= E[log-odds]
%       x((2^K)+1:2^(K+1))= log V[log-odds]
%   - phi: phi(1) = (log-) temperature and phi(2) = bias
%   - u: u(1:K)= sequence of past outcomes
%   - in: depth of sequence learning
% OUT:
%   - gx: P(y_t=1|y_{t-1})

a = 0.36; % for E[s(x)] when x~n(mu,Sig)
K = in.K; % sequence depth
yb = u(2:K+1); % previous outcomes

if ~any(isnan(u)) % trial OK
    if K >0
        indSeq = bin2dec(num2str(yb'))+1; % index of sequence of previous outcomes
    else
        indSeq = 1;
    end
    m = x(indSeq);
    v = exp(x((2^K)+indSeq));
    gx = sig(phi(2)+exp(phi(1)).*m./sqrt(1+a*v)); % E[sigm(log-odds of P(y))]
else
    gx = 0.5;
end

 