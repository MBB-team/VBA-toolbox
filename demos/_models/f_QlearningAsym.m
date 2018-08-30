function  [fx, dfdx, dfdp] = f_QlearningAsym (x, P, u, in)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [fx, dfdx, dfdp] = f_QlearningAsym (x, P, u, in)
% Reinforcement-learning evolution function for a n-armed bandit task with
% asymmetric learning for positive and negative feedbakcs
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
%        (3:4) unused
%	- in: [useless]
% OUT:
%   - fx: updated action values
%   - dfdx/dfdP: gradients for VBA inversion

% /////////////////////////////////////////////////////////////////////////

n = numel(x);

% In the case there was no feedback, do nothing
% =========================================================================

if isnan(u(2))
    fx = x;
    dfdx = eye(n);
    dfdp = zeros(numel(P),n);
    return;
end

% Apply delta-rule to update action values
% =========================================================================

% get experimental conditions
cueIdx = u(3 : 4);
prevActionIdx = cueIdx(u(1) + 1); % action 0 is first index
feedback = u(2);
    
% start with previous values
fx = x; 

% prediction error
delta = feedback - x(prevActionIdx);

% asymmetric learning rate 
alpha = VBA_sigmoid(P(1) + sign (delta) * P(2)); % [-Inf,Inf] -> [0 1]

% update Q-value
fx(prevActionIdx) = x(prevActionIdx) + alpha*delta; % update chosen value


% Compute evolution function's gradient
% =========================================================================
% This is not necessary, as the toolbox will approximate those gradients if
% needed. However, providing the analytical gradient can higly speed-up the
% inversion.


% derivative w.r.t hidden state
% -------------------------------------------------------------------------
dfdx = eye(n);
dfdx(prevActionIdx,prevActionIdx) = 1 - alpha;

% derivative w.r.t parameters
% -------------------------------------------------------------------------
dfdp = zeros(2,n);
dfdp(1,prevActionIdx) = alpha * (1 - alpha) * delta;
dfdp(2,prevActionIdx) = alpha * (1 - alpha) * sign(delta) * delta;
