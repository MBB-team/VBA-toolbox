function [priors] = VBA_defaultMFXPriors (dim)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [priors] = VBA_defaultMFXPriors (dim)
% returns a default priors structure for hiearchical (random effect)
% inference
%
% This function sets the default parameters associated to the prior pdf
% of the group level parameters -- that group mean and variance of Phi, 
% Theta, and X0 -- which are required for the inversion.
% For the population mean the default prior is the unit multivariate 
% Gaussian N(0,Id).
% For the population variance the prior over precision is the Gamma 
% distribution Ga(1,1).
%
% IN:
%   - dim: a structure containing the dimensions of the model variables, ie.
%       .n_phi: number of observation parameters
%       .n_theta: number of evolution parameters
%       .n: number of hidden states
%
% OUT:
%   - priors: a structure containing the parameters of the prior pdf, ie:
%       .muPhi, .muTheta, .muXO: prior mean on population average of the
%       observation (1 x n_phi) and evolution (1 x nTheta) parameters, and 
%       inital state (1 x n).
%       .SigmaPhi, .SigmaTheta, .SigmaX0: prior variance on population
%       average of the observation (n_phi x n_phi) and evolution 
%       (nTheta x nTheta) parameters, and inital state (n x n).
%       .a_vPhi, .a_vTheta, .a_vX0: prior shape param on population variance
%       .b_vPhi, .b_vTheta, .b_vX0: prior scale param on population variance
%
% /////////////////////////////////////////////////////////////////////////

priors = struct ();

% model parameters
% -------------------------------------------------------------------------
% pdf of the observation parameters
priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = eye(dim.n_phi);
priors.a_vPhi = ones(dim.n_phi,1);
priors.b_vPhi = ones(dim.n_phi,1);

% pdf of the evolution parameters
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaTheta = eye(dim.n_theta);
priors.a_vTheta = ones(dim.n_theta,1);
priors.b_vTheta = ones(dim.n_theta,1);

% initial hidden states
% -------------------------------------------------------------------------
% pdf for the initial state
priors.muX0 = zeros(dim.n,1);
priors.SigmaX0 = eye(dim.n);
priors.a_vX0 = ones(dim.n,1);
priors.b_vX0 = ones(dim.n,1);

