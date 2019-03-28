function  [gx] = g_Kalman4bandits (x, P, ~, in)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [gx] = g_Kalman4bandits (x, P, u, in)
% softmax decision rule for N-armed bandit task as in Daw et al. 2006
%
% /////////////////////////////////////////////////////////////////////////

% Get parameter values
% -------------------------------------------------------------------------
% inverse temperature
beta = exp (P(1)); % exp: [-Inf,Inf] -> [0 Inf]
% bonus to exploration
phi = P(2); 

% Extract expectation and variance of each bandit
% -------------------------------------------------------------------------
nBandits = in.nBandits;
mu = x(1 : nBandits); 
sigma2 = x(nBandits + (1 : nBandits)) .^ 2;

% apply softmax
% -------------------------------------------------------------------------
gx = VBA_softmax(beta * (mu + phi * sigma2));