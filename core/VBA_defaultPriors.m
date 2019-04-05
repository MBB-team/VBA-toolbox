function [priors] = VBA_defaultPriors (dim, options)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [priors] = VBA_defaultPriors (dim, options)
% returns a default priors structure
%
% This function sets the default parameters associated to the prior pdf
% which are required for the model inversion.
% For evolution (theta), observation (phi) and initial state (X0), the
% default prior is the unit multivariate Gaussian N(0,Id).
% For the measurement noise (sigma) and the state stochastic innovation
% (alpha), the (Jeffrey) prior over precision is the Gamma distribution
% Ga(1,1).
%
% IN:
%   - dim: a structure containing the dimensions of the model variables, ie.
%       .n_phi: number of observation parameters
%       .n_theta: number of evolution parameters
%       .n: number of hidden states
%       .n_t number of timepoints
%   - options: a structure with the field
%       .source: structure defining the type of observations
%
% OUT:
%   - priors: a structure containing the parameters of the prior pdf, ie:
%       .muPhi: a n_phix1 vector containing the prior mean of Phi, the
%       observation parameters
%       .muTheta: a n_thetax1 vector containing the prior mean of Theta,
%       the evolution parameters
%       .muX0: a nx1 vector containing the prior mean of the hidden-states
%       initial condition
%       .SigmaPhi: n_phixn_phi prior covariance matrix of Phi
%       .SigmaTheta: n_thetaxn_theta prior covariance matrix of Theta
%       .SigmaX0: nxn prior covariance matrix of X0
%       .a_sigma / .b_sigma: the shape and scale parameters of the prior
%       Gamma pdf upon the measurement noise precision
%       .a_alpha / .b_alpha: the shape and scale parameters of the prior
%       Gamma pdf upon the stochastic innovations precision
%
% /////////////////////////////////////////////////////////////////////////


% check for sources for backward compatibility
if ~isfield(options,'sources')
    error('*** VBA_defaultPriors: please specify the observation distribution type using options.sources.');
end

% model parameters
% -------------------------------------------------------------------------
% prior Gaussian pdf of the observation parameters
priors.muPhi = zeros(dim.n_phi, 1);
priors.SigmaPhi = eye(dim.n_phi);

% prior Gaussian pdf of the evolution parameters
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaTheta = eye(dim.n_theta);

% measurement noise precision, for Gaussian sources
% -------------------------------------------------------------------------
gsi = find([options.sources.type] == 0);
ngs = numel(gsi);

% prior Gamma pdf (Jeffrey)
priors.a_sigma = 1e0 * ones(1, ngs);
priors.b_sigma = 1e0 * ones(1, ngs);

% Covariance structure
priors.iQy = cell (dim.n_t, length (gsi));
for i = 1 : length (gsi)
    nY = length (options.sources(gsi(i)).out);
    priors.iQy(:, i) = {eye(nY)};
end

% hidden states
% -------------------------------------------------------------------------
% prior Gaussian pdf for the initial state
priors.muX0 = zeros(dim.n,1);
priors.SigmaX0 = eye(dim.n);
    
% state noise precision
if dim.n > 0
    % singular prior Gamma pdf (deterministic system)
    priors.a_alpha = Inf;
    priors.b_alpha = 0;
    % Covariance structure
    priors.iQx = cell(dim.n_t,1);
    priors.iQx(:) = {speye(dim.n)};
else
    priors.a_alpha = [];
    priors.b_alpha = [];
    priors.iQx = [];
end


