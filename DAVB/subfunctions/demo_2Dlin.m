% Demo for partially observable oscillatory system.
% This demo inverts a model of a linear oscillatory system, which is
% observed through a nonlinear sigmoid function.
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

clear variables
close all

% Choose basic settings for simulations
n_t = 5e2;
delta_t = 1e-1;
f_fname = @f_lin2D;
g_fname = @g_Id;

u       = [];


% Parameters of the simulation
alpha   = 1e2;
sigma   = 1e4;
theta   = 1;
phi     = [];

% Build priors for model inversion
priors.muX0 = [-2;-2];
priors.SigmaX0 = 1e-0*eye(2);
priors.muTheta = 0*ones(1,1);
priors.SigmaTheta = 1e-1*eye(1);
priors.a_alpha = 1e2;
priors.b_alpha = 1e1;
priors.a_sigma = 1e4;
priors.b_sigma = 1e1;

% Build options and dim structures for model inversion
options.priors      = priors;
options.inF.deltat = delta_t;
options.backwardLag  = 2;
dim.n_theta         = 1;
dim.n_phi           = 0;
dim.n               = 2;

% options.checkGrads = 1;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
disp('--paused--')
pause

dbstop if error

% Call inversion routine
% [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options);
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = comparePredictions(...
        n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end


