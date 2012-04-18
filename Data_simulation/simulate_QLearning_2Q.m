function [A,r,Q,p] = simulate_QLearning_2Q(Q0,alpha,beta,V)
% This function simulates a model learning 2 Qvalues through a Rescorla
% Wagner update rule and choosing one them using a softmax decision rule.
% 
% IN :
%   - Q0: 2x1 initial Qvalues
%   - alpha: learning rate
%   - beta: inverse temperature of softmax decision rule
%   - V: 2*N the values of rewards for each alternative (computed in advance)
% OUT:
%   - A: 1xN the sequence of actions of the model (binary)
%   - r: 1xN the obtained reward
%   - Q: 2xN the vector of Qvalues
%   - p: 2xN the probability of selection of each alternative
%------------------------------------------------------------------------
% Vincent Adam 02/04/2012
%------------------------------------------------------------------------


N = length(V); % Number of trials
A = zeros(1,N); % chosen action time series
p = zeros(2,N); % action likelihood
r = zeros(1,N); % obtained feedback
Q = zeros(2,N+1); % Q-values

% First trial
% action emission
p1 = 1/(1+exp(beta*(Q0(2)-Q0(1)))); % Probability of choosing choice 2
p(:,1) = [p1;1-p1];
a = sampleFromArbitraryP(p(:,1)',[0,1],1);   % Sample action
A(1)=a;
r(1) = V(a+1,1); % feedback
Q(:,1)= Q0; % Q-values update
Q(a+1,1)=Q0(a+1)+alpha*(r(1)-Q0(a+1)); % Learning from feedback


% Simulating data
for t = 2 : N
    % action emission
    p1 = 1/(1+exp(beta*(Q(2,t-1)-Q(1,t-1)))); % Probability of choosing choice 2
    p(:,t) = [p1;1-p1];
    a = ~sampleFromArbitraryP(p(:,t)',[0,1],1);   % Sample action
    A(t)=a;
    r(t) = V(a+1,t);    % feedback
    Q(:,t)= Q(:,t-1);    % Q-values update
    Q(a+1,t)=Q(a+1,t-1)+alpha*(r(t)-Q(a+1,t-1)); % Learning from feedback
    
end
Q = Q(:,1:end-1); % Deleting last update not used for action selection