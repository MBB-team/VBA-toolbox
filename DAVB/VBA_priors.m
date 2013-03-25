function [priors] = VBA_priors(dim,options)
% checks (and fills in default) priors for the sDCM generative model
% function [priors] = VBA_priors(dim)
% This function sets the default parameters associated to the prior pdf
% which are required by the NL state-space model.
% IN:
%   - dim: a structure containing the dimensions of the model variables
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


if isfield(options,'extended') && options.extended
    [priors] = VBA_priors_extended(dim,options);
    return
end

% prior Gamma pdf of the measurement noise (Jeffrey)
if ~options.binomial
    priors.a_sigma = 1e0;
    priors.b_sigma = 1e0;
end

% Covariance structure: measurement noise precision matrices
priors.iQy = cell(dim.n_t,1);
for t=1:dim.n_t
    priors.iQy{t,1} = eye(dim.p);
end

% prior Gaussian pdf of the observation parameters
if dim.n_phi > 0
    priors.muPhi = zeros(dim.n_phi,1);
    priors.SigmaPhi = eye(dim.n_phi);
else
    priors.muPhi = [];
    priors.SigmaPhi = [];
end

% prior Gaussian pdf of the evolution parameters
if dim.n_theta > 0
    priors.muTheta = zeros(dim.n_theta,1);
    priors.SigmaTheta = eye(dim.n_theta);
else
    priors.muTheta = [];
    priors.SigmaTheta = [];
end

if dim.n > 0
    % prior Gaussian pdf of the hidden states
    priors.muX0 = zeros(dim.n,1);
    priors.SigmaX0 = eye(dim.n);
    % prior Gamma pdf of the state noise (Jeffrey)
    priors.a_alpha = 1e0;
    priors.b_alpha = 1e0;
    % Covariance structure: state noise precision matrices
    priors.iQx = cell(dim.n_t,1);
    for t=1:dim.n_t
        priors.iQx{t} = eye(dim.n);
    end
else
    priors.muX0 = [];
    priors.SigmaX0 = [];
    priors.a_alpha = [];
    priors.b_alpha = [];
    priors.iQx = [];
end


