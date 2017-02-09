function [fx]=evolution0bisND(x,theta, u,in)
% 0-ToM's evolution function (without doubled hidden states)
% function [fx]=evolution0bisND(x,theta, u,in)
% 0-ToM is simply tracking the log-odds of P(o=1), where o is the
% opponent's action. This variable is updated according to a Laplace-Kalman
% filter, yielding 2 sufficient statistics, mu and Sig. In this scheme, the
% only evolution param (theta) is 0-ToM's prior volatity about her
% opponent's log-odds. 
% IN:
%   - x: sufficient statistics of log-odds of P(o=1):
%       x(1)= E[log-odds]
%       x(2)= log V[log-odds] (log-scale for numerical reasons)
%   - theta: 0-ToM's prior volatity
%   - u: last opponent's move (o)
%   - in: [useless]
% OUT:
%   - fx: updated sufficient statistics of log-odds of P(o=1)

% theta is also transform through 1.4* sigmoid() so that no explosion is
% created when simulated by higher ToM-levels

if isempty(u)||isnan(u(1)) % missed trial
    fx=x; % no update
else    
    mu=x(1); % current E[log-odds]
    Sig=exp(x(2)); % current V[log-odds]
    % Laplace-Kalman update rule
    volatility = 1.4*sigmoid(theta(1)); % for numerical reasons (when simulated by k-ToM with k>0)
    SigT=1/((1/(volatility+Sig))+sigmoid(mu)*(1-sigmoid(mu)));
%     SigT=1/(1/(exp(theta)+Sig)+sigmoid(mu)*(1-sigmoid(mu))+1e-2);
    muT=mu+SigT*(u(1)-sigmoid(mu));
    % wrap-up
    fx=[muT;log(SigT)];
end    