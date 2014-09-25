% Demo for stochastic system with binomial output
% This demo inverts a model of an AR system, which is observed through a
% nonlinear sigmoid function.

clear variables
close all

% Choose basic settings for simulations
n_t = 1e2;
f_fname = @f_AR;
g_fname = @g_probit;

u       = [];


% Parameters of the simulation
n = 1; % # hidden states
alpha   = 1e1;
theta   = [];
phi     = 0;
x0 = zeros(n,1);

% Build priors for model inversion
priors.muX0 = zeros(n,1);
priors.SigmaX0 = 1e0*eye(n);
priors.muPhi = 0*ones(1,1);
priors.SigmaPhi = 1e0*eye(1);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;


% Build options and dim structures for model inversion
options.priors      = priors;
options.inG.bias   = randn(4,1);
options.backwardLag  = 4;
options.binomial = 1;
dim.n_theta         = 0;
dim.n_phi           = 1;
dim.n               = n;

% options.checkGrads = 1;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,[],options,x0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)

% Call inversion routine
% [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options);
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,[])

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = comparePredictions(...
        n_t,theta,phi,u,alpha,[],options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end


