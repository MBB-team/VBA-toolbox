function [posterior, out] = demo_logisticRegression()
% // VBA toolbox //////////////////////////////////////////////////////////
% [posterior, out] = demo_logisticRegression()
%
% this function exemplifies the use of the logistic model to perform a
% regression on binary data

% Specify the logistic model
% =========================================================================

% number of observations
n = 128;

% number of predictors
nRegressor = 4;

% design matrix (each regressor is entered as a line)
X = [ ones(1,n)             ;   % constant (intercept)
      randn(nRegressor-1,n) ] ; % random predictors

% Prepare options structure
% =========================================================================

% model definition
% -------------------------------------------------------------------------
g_fname = @g_logistic  ;
u = X ;        

% define model dimensions
% -------------------------------------------------------------------------
dim.n_phi   = nRegressor;  % nb of observation parameters
dim.n_theta = 0;         % nb of evolution parameters
dim.n       = 0;                 % nb of
options.dim = dim;

% Specify distribution for binary data (default is gaussian)
% -------------------------------------------------------------------------
options.sources.type = 1; 

% Simulate artificial data
% =========================================================================

% Choose parameters for simulations (weight of predictors)
phi = ones(nRegressor,1); % basic example

% Simulated observations 
[y,x,x0,eta,e] = VBA_simulate (n,[],g_fname,[],phi,u,[],[],options);

% display time series of hidden states and observations
figure,
plot(y,'ro') % simulated data
hold on;
plot((y-e))  % model prediction = data - residual

% Run the estimation  
% =========================================================================

% Call inversion routines
[posterior,out] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);

% Display results
displayResults(posterior,out,y,[],[],[],phi,[],[]);