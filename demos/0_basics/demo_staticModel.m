function [posterior, out] = demo_staticModel ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_staticModel ()
% Demo of simple linear Bayesian regression
%
% The problem is to fit observations with a linear combination
% of regressors. Formally, this can be written as:
%
% $$ y = X.\beta + \vareps $$
%
% where $y$ is the set of observations, $X$ is the design matrix containing
% the regressors, $\beta$ is a vector of parameters defining the respective 
% weight of the regressors, and $\vareps ~ N(0,\sigma^2)$ is some measurment
% noise.
% The objective is thus to estimate $\beta$ and $\sigma^2$ from the data
% given X.
%
% /////////////////////////////////////////////////////////////////////////

%% Specify the linear model
% =========================================================================

% number of observations
N = 100;

% regressors
k = 1; % this is a simple regression. Set k > 1 for a multiple regression model
regressors = randn(k, N);

% observation function that links inputs (regressors), parameters (beta
% weights) and observations
function [gx] = g_linearRegression (~, betas, regressors, ~)
    % predictied observations
    gx = betas' * regressors ;
end

%% Simulate data
% =========================================================================
% 'true' weights of the regressors
betas = rand(k,1) ; 

% 'true' observation noise
sigma2 = 5;

% simulate data
% this applies g_linearRegression to all columns of the inputs (regressors) 
% with the given parameters beta and sigma
[y] = VBA_simulate (N,[],@g_linearRegression,[],betas,regressors,Inf,sigma2,struct,[]);


%% Inversion
% =========================================================================
% specify dimensions of the problem
% + number of observation parameters, ie. here nuber of betas
dim.n_phi = k;
% number of inputs, ie. here number of regressors
dim.u = k;
% + size of observations (number of predictions per input)
dim.p = 1;

% call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,regressors,[],@g_linearRegression,dim,struct);

%% Results
% =========================================================================
% estimated beta weights
fprintf('Estimated beta weights: ');
fprintf('%3.2f ', posterior.muPhi);
fprintf('\n');

% fit
fprintf('Percentage of explained variance (R2): %03.2f\n', out.fit.R2);

% plot results vs. true parameter values
displayResults(posterior,out,y,[],[],[],betas,[],sigma2);

end

