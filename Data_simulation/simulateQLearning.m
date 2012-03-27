function [A,r,Q,p] = simulateQLearning(Q0,alpha,beta,V)

N = length(V);
A = zeros(1,N); % chosen action time series
p = zeros(2,N); % action likelihood
r = zeros(1,N); % obtained feedback
Q = zeros(2,N+1); % Q-values

% First trial
% action emission
p2 = 1/(1+exp(beta*(Q0(2)-Q0(1)))); % Probability of choosing choice 2
p(:,1) = [1-p2;p2];
a = sampleFromArbitraryP(p(:,1)',[1,0],1);   % Sample action
A(1)=a;
% feedback
r(1) = V(a+1,1);
% Q-values update
Q(:,1)= Q0;

Q(a+1,1)=Q0(a+1)+alpha*(r(1)-Q0(a+1)); % Learning from feedback


% Simulating data
for t = 2 : N
    
    % action emission
    p2 = 1/(1+exp(beta*(Q(2,t-1)-Q(1,t-1)))); % Probability of choosing choice 2
    p(:,t) = [1-p2;p2];
    a = sampleFromArbitraryP(p(:,t)',[1,0],1);   % Sample action
    A(t)=a;
    % feedback
    r(t) = V(a+1,t);
    % Q-values update
    Q(:,t)= Q(:,t-1);
    
    Q(a+1,t)=Q(a+1,t-1)+alpha*(r(t)-Q(a+1,t-1)); % Learning from feedback
    
end

Q = Q(:,1:end-1); % Deleting last update not used for action selection