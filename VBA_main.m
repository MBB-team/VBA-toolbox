function [posterior,out] = VBA_main(y,M)
% inverts nonlinear state-space model M given data y
% function [posterior,out] = VBA_main(y,M)
% IN:
%   - y: observed data matrix
%   - M: model structure, with fields:
%       .u: input to the system {[]}
%       .f_fname: handle of the evolution mapping
%       .g_fname: handle of the observation mapping
%       .dim: structure of model variables dimensions
%       .options: structure with optional i/o
% OUT:
%   - posterior: a structure variable whose fields contains the sufficient
%   statistics (typically first and second order moments) of the
%   variational approximations of the posterior pdfs over the
%   observation/evolution/precision parameters and hidden-states time
%   series. Its fields are:
%   	.muX: posterior mean of the hidden states X (nxn_t matrix)
%       .SigmaX: covariance matrices of the variational posterior pdf of
%       the dynamic hidden-states. Using the Kalman-Rauch marginalization
%       procedure, this is further divided into:
%           SigmaX.current{t} : the instantaneous covariance matrix at t
%           SigmaX.inter{t} : the lagged covariance matrix between instants
%           t and t+1
%       .muX0: posterior mean of the hidden-states initial condition, ie
%       before the first observation (nx1 vector)
%       .SigmaX0: covariance matrix of the Gaussian df over the
%       hidden-states initial condition (nxn matrix)
%   	.muTheta: posterior mean of the evolution parameters (vector)
%       .SigmaTheta: covariance matrix of the variational posterior pdf of
%       the static evolution parameters
%   	.muPhi: posterior mean of the observation parameters (vector)
%       .SigmaPhi: covariance matrix of the variational posterior pdf of
%       the static observation parameters
%       .a_alpha / .b_alpha: shape and scale parameters of the variational
%       posterior pdf of the stochastic innovations precision
%       .a_sigma / .b_sigma: shape and scale parameters of the variational
%       posterior pdf of the measurement noise precision
%   - out: a structure variable containing the fields...
%       .CV: convergence flag (0 if the algorithm has stopped because it
%       reached the options.MaxIter termination condition)
%       .F: the free energy associated with the inversion of the model
%       .M: the model structure with all fields filled in for book keeping
%       .it: the number of iterations which have been required for reaching
%       the convergence criteria
%       .suffStat: a structure containing internal variables that act as
%       sufficient statistics for the VB updates (e.g. predicted data ...)

if ~isfield(M,'options')
    M.options = [];
end
if ~isfield(M,'u')
    M.u = [];
end

[posterior,out] = VBA_NLStateSpaceModel(y,M.u,M.f_fname,M.g_fname,M.dim,M.options);