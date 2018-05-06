function  [gx, dgdx, dgdp] = g_QLearning (x, P, u, in)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [gx, dgdx, dgdP] = g_QLearning (x, P, u, in)
% softmax decision rule for Q-learning (N-armed bandit task)
%
% IN:
%   - x: Q-values
%   - P: (1) inverse (log-) temperature 
%        (2) bias
%   - u: (idx(1)) index of 0 coded cue
%        (idx(2)) index of 1 coded cue
%   - in: 
%       - idx: position of inputs indicating Q-values to use
% OUT:
%   - gx : P(y = 1|x)

% /////////////////////////////////////////////////////////////////////////

% Get parameter values
% =========================================================================

% inverse temperature
beta = exp (P(1)); % exp: [-Inf,Inf] -> [0 Inf]

% offset
try
    const = P(2);
catch 
    const = 0;
end

% Behavioural prediction
% =========================================================================

% get idx of Q-value to use
if numel (x) == 2
    % if only two arms, directly takes the two Q-values
    idx = [1 2]; 
else
    % otherwise, by default the last 2 inputs should provide the curent
    % cues
    idx = u(end-1:end);
end

% Compute Value differential
dQ = x(idx(2))-x(idx(1));

% make prediction
gx = VBA_sigmoid(beta * dQ + const);


% Compute gradients
% =========================================================================

% w.r.t hidden state
% -------------------------------------------------------------------------
dgdx = zeros (numel (x), 1);
dgdx(idx) = [-1 1] * beta * gx * (1 - gx);


% w.r.t parameters
% -------------------------------------------------------------------------
dgdp = zeros (numel (P), 1);

% temperature
dgdp(1) = beta * dQ * gx * (1 - gx);

% offset
if numel (P) > 1
    dgdp(2) = gx * (1 - gx);
end


