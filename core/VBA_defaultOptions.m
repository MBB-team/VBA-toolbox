function [options] = VBA_defaultOptions ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [options] = VBA_defaultOptions ()
% returns a default options structure
%
% OUT:
%   - options: a structure containing the parameters of the inversion:
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

options = struct(...
    'decim'      , 1     , ...     % Micro-time resolution
    'microU'     , 0     , ...     % Micro-resolution input
    'inF'        , []    , ...     % Optional (internal) parameters of the evolution function
    'inG'        , []    , ...     % Optional (internal) parameters of the observation function
    'checkGrads' , 0     , ...     % Check analytical gradients against numerical gradients
    'updateX0'   , 1     , ...     % Flag for initial conditions update
    'updateHP'   , 1     , ...     % Flag for hyperparameters update
    'initHP'     , 1     , ...     % Flag for hyperparameters initialization
    'backwardLag', 1     , ...     % time lag of the short-sighted backward-pass
    'MaxIter'    , 32    , ...     % Maximum number of iterations
    'MaxIterInit', 32    , ...     % Maximum number of iterations for the initialization
    'MinIter'    , 0     , ...     % Minimum number of iterations
    'TolFun'     , 2e-2  , ...     % Minimum change in the free energy
    'DisplayWin' , 1     , ...     % VB display window
    'gradF'      , 0     , ...     % Gauss-Newton ascent on free (1) or variational (0) energy
    'GnMaxIter'  , 32    , ...     % Gauss-Newton coordinate ascent maximum number of iterations:
    'GnTolFun'   , 1e-5  , ...     % Gauss-Newton coordinate ascent threshold on relative variational energy:
    'GnFigs'     , 0     , ...     % Gauss-Newton inner-loops figures
    'verbose'    , 1     , ...     % matlab window messages
    'OnLine'     , 0     , ...     % On-line version (true when called from VBA_OnLineWrapper.m)
    'delays'     , []    , ...     % delays
    'kernelSize' , 16    , ...     % max lag of Volterra kernels
    'nmog'       , 1     , ...     % split-Laplace VB?
    'UNL'        , 0     , ...     % un-normalized likelihood?
    'UNL_width'  , 4     , ...     % for partition function estimation
    'UNL_ng'     , 64      ...     % for partition function estimation
);