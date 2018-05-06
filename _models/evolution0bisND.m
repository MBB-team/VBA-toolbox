function [fx] = evolution0bisND(x,theta,u,in)
% 0-ToM's evolution function (without doubled hidden states)
% function [fx] = evolution0bisND(x,theta,u,in)
% 0-ToM is simply tracking the log-odds of P(o=1), where o is the
% opponent's action. This variable is updated according to a Laplace-Kalman
% filter, yielding 2 sufficient statistics, m and V. In this scheme, the
% only evolution param (theta) is 0-ToM's prior volatity about her
% opponent's log-odds. 
% IN:
%   - x: sufficient statistics of log-odds of P(o=1):
%       x(1)= E[log-odds]
%       x(2)= log V[log-odds] (log-scale for numerical reasons)
%   - theta: 0-ToM's prior (log-) volatity
%   - u: u(1)= last opponent's move (o)
%   - in: [useless here]
% OUT:
%   - fx: updated sufficient statistics of log-odds of P(o=1)

if isempty(u)||isnan(u(1)) % missed trial
    fx = x; % no update
else % trial OK
    % -- deal with superceding competition with other generative models --
    % [See, e.g., f_metaTom.m]
    try
        w = inF.metaweight;
    catch
        w = 1;
    end
    % -- learning rule --
    m0 = x(1); % current E[log-odds]
    V0 = exp(x(2)); % current V[log-odds]
    p0 = VBA_sigmoid(m0); % current estimate of P(o=1)
    volatility = exp(theta(1));
    V = 1./((1./(volatility+V0))+w*p0*(1-p0)); % updated V[log-odds]
    m = m0 + w*V*(u(1)-p0); % updated E[log-odds] (Laplace-Kalman update rule)
    % wrap-up
    fx = [VBA_sigmoid(VBA_sigmoid(m),'inverse',true);log(V)]; % for numerical purposes
end    