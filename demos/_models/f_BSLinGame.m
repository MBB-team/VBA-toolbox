function [ fx ] = f_BSLinGame(x,theta,u,in)
% evolution function for BSL learner engaging in dyadic games
% [ fx ] = f_BSLinGame(x,theta,u,in)
% BSL is simply tracking the log-odds of P(o_t=1|o_{t-1}), where o is the
% opponent's binary action. This variable is updated according to a
% Laplace-Kalman filter, yielding 2 sufficient statistics (m and V) per
% combination of past outcome. BSL can learn sequences of arbitrary depth
% (K). For example, if K=1, then BSL tracks 2 probabilities, namely: 
% P(o_t=1|o_{t-1}=1) and P(o_t=1|o_{t-1}=0). More generally, BSL tracks 2^K
% probabilities. In this scheme, the only evolution param (theta) is BSL's
% prior volatity about the log-odds.
% Note: unsampled sequences will eventually be "forgotten", since the
% prediction step in the Laplace-Kalman update will dilute any previously
% sampled evidence.
% IN:
%   - x: sufficient statistics of log-odds of P(o=1):
%       x(1:2^K)= E[log-odds]
%       x((2^K)+1:2^(K+1))= log V[log-odds]
%   - theta: BSL's prior volatity
%   - u: sequence of past outcomes:
%       u(1)= opponent's last move
%       u(2)= learner's last move
%       u(3:K+2) = sequence of K past opponent's moves
%   - in: depth of sequence learning
% OUT:
%   - fx: updated sufficient statistics of log-odds of P(o=1)

if VBA_isWeird (u) % e.g., 1st trial
    fx = x;
    return
end

u(2) = []; % remove agent's previous move
[fx] = f_BSL(x,theta,u,in);