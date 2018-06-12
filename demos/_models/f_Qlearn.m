function  [fx,dfdx,dfdP] = f_Qlearn(x,P,u,in)
% Reinforcement-learning model for a 2-armed bandit task
% function  [fx,dfdx,dfdP] = f_Qlearn(x,P,u,in)
% An RL agent learns by trial and error. A bandit task is such that, after
% each action, the agent receives a feedback (reward if positive,
% punishment if negative). The RL agent updates its action values as
% follows:
% V(chosen action) = V(chosen action) + alpha*(feedback-V(chosen action))
% V(unchosen action) = V(unchosen action)
% IN: 
%	- x: action values
%	- P: (invsigmoid-) learning rate
%	- u: previous action (u(1)) and feedback (u(2))
%	- in: [useless]
% OUT:
%   - fx: updated action values
%   - dfdx/dfdP: gradients for VBA inversion

alpha = VBA_sigmoid(P(1)); % learning rate
a = 2-u(1); % index of agent's last chosen action
r = u(2); % feedback
fx = x; % identity mapping
fx(a) = x(a) + alpha*(r-x(a)); % update chosen value

% compute gradients
if a == 1
    dfdx = [1-alpha, 0;
            0, 1];
    dfdP = [alpha*(1-alpha)*(r-x(a)),0];
elseif a == 2                   
    dfdx = [1, 0;
            0, 1-alpha];
    dfdP = [0,alpha*(1-alpha)*(r-x(a))];
end

% dfdx = dfdx';
% dfdP = dfdP';