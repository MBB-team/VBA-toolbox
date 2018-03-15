function [posterior, out] = demo_GLM()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_GLM()
% Demo of multiple regression (General Linear Model) inference
%
% The problem is to fit observations with a linear combination
% of regressors. Formally, this can be written as:
%
% $$ y = X.\beta + \vareps $$
%
% where $y$ is the set of observations, $X$ is the design matrix containing
% the regressors, $\beta$ is a vector of parameters defining the respective 
% weight of the regressors, and $\vareps ~ N(0,\sigma^)$ is some measurment
% noise.
% The objective is thus to estimate $\beta$ and $\sigma^2$ from the data
% given X.
%
% /////////////////////////////////////////////////////////////////////////

% Specify the linear model
% =========================================================================

% number of observations
n = 8;

% design matrix (each regressor is entered as a column)
X = [ ones(n,1)    ... % constant (intercept)
      (1:n)'       ] ; % linear trend

% Generate artificial data
% =========================================================================

% Draw random effect strength
nRegressor = size(X,2);
b = randn(nRegressor,1);

% Simulate random data with precision sigma
sigma = 1;    % precision = 1 / variance
y0    = X*b;  % model prediction
y = y0 + randn(n,1)/sigma;

% Prepare options structure for the inversion
% =========================================================================

% model definition
% -------------------------------------------------------------------------
g_fname = @g_GLM  ;
options.inG.X = X ; % will be used by g_GLM

% define priors
% -------------------------------------------------------------------------
% % default priors are N(0,1) for all parameters. 
% % uncomment the following lines to change the default
% % prior mean on observation params (beta values)
%   priors.muPhi = zeros(nRegressor,1);    
% % prior covariance on observation params
%   priors.SigmaPhi = 1e0*eye(nRegressor); 
% % hyperpriors
%   priors.a_sigma = 1; % Jeffrey's priors
%   priors.b_sigma = 1; % Jeffrey's priors
% % include priors in options structure
%   options.priors = priors;        

% define model dimensions
% -------------------------------------------------------------------------
dim.n_phi = nRegressor;  % nb of observation parameters
dim.n_theta = 0;         % nb of evolution parameters
dim.n=0;                 % nb of hidden states 

% Run the estimation 
% =========================================================================

% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);

% Display results
displayResults(posterior,out,y,[],[],[],b,[],sigma);
