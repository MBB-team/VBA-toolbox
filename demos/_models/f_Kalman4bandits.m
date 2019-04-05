function  [fx] = f_Kalman4bandits (x, P, u, in)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [fx] = f_Kalman4bandits (x, P, u, in)
% 
% Kalman filter learning rule for the N-armed bandit task as in 
% Daw et al. 2006
%
% /////////////////////////////////////////////////////////////////////////

% rename entries
% =========================================================================
% task dimension
nBandits = in.nBandits;

% inputs
choice = u(1 : nBandits);
feedback = u(nBandits + 1);

% parameters
lambda = P(1);
theta = P(2);
sigma_d = P(3);
sigma_o = P(4);

% states
mu = x(1 : nBandits);
sigma2 = x(nBandits + (1 : nBandits)) .^ 2;

% In the case there was no feedback (first trial), do nothing
% =========================================================================
if isnan (feedback)
    fx = x;
    return;
end

% Apply delta-rule to update action values
% =========================================================================
% chosen option
% -------------------------------------------------------------------------
% index
chosen = choice == 1;

% prediction error
delta = feedback - mu(chosen);

% update
kappa = sigma2(chosen) / (sigma2(chosen) + sigma_o^2);
mu(chosen) = mu(chosen) + kappa * delta;
sigma2(chosen) = (1 - kappa) * sigma2(chosen);

% unchosen options
% -------------------------------------------------------------------------
% index
unchosen = choice == 0;

% update
sigma2(unchosen) = lambda^2 * sigma2(unchosen) + sigma_d^2;
mu(unchosen) = lambda * mu(unchosen) + (1 - lambda) * theta;

% store
% =========================================================================
fx = [mu; sqrt(sigma2)];
