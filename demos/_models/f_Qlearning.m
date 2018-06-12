function  [fx,dfdx,dfdp] = f_Qlearning(x,P,u,in)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [fx,dfdx,dfdP] = f_Qlearning(x,P,u,in)
% Reinforcement-learning evolution function for a n-armed bandit task
%
% An RL agent learns by trial and error. A bandit task is such that, after
% each action, the agent receives a feedback (reward if positive,
% punishment if negative). The RL agent updates its action values as
% follows:
% V(chosen action) = V(chosen action) + alpha*(feedback-V(chosen action))
% V(unchosen action) = V(unchosen action)
% IN: 
%	- x: action values (n x 1)
%	- P: learning rate (will be sigmoid transformed)
%	- u: (1) previous action 
%        (2) feedback received for previous action
%	- in: [useless]
% OUT:
%   - fx: updated action values
%   - dfdx/dfdP: gradients for VBA inversion

% /////////////////////////////////////////////////////////////////////////


% Get parameter values
% =========================================================================
% Some of the model parameters can only take values within a certain range.
% However, the VBA toolbox can only perform estimation of unbounded values.
% To solve this, we transform the values passed by VBA to map them on the
% acceptable range for each parameter.

% learning rate 
alpha = VBA_sigmoid(P); % [-Inf,Inf] -> [0 1]


% Apply delta-rule to update action values
% =========================================================================

% get experimental conditions
prevActionIdx = u(1)+1; % action 0 is first index
feedback = u(2);

% start with previous values
fx = x; 

% udpdate previous action value
delta = feedback - x(prevActionIdx);
fx(prevActionIdx) = x(prevActionIdx) + alpha*delta; % update chosen value


% Compute evolution function's gradient
% =========================================================================
% This is not necessary, as the toolbox will approximate those gradients if
% needed. However, providing the analytical gradient can higly speed-up the
% inversion.

n = numel(x);

% derivative w.r.t hidden state
% -------------------------------------------------------------------------
dfdx = eye(n);
dfdx(prevActionIdx,prevActionIdx) = 1 - alpha;

% derivative w.r.t parameters
% -------------------------------------------------------------------------
dfdp = zeros(1,n);
dfdp(prevActionIdx) = alpha*(1-alpha)*delta;
