% Demo for partially observable oscillatory system.
% This demo inverts a model of a linear oscillatory system, which is
% observed through a nonlinear sigmoid function.

clear variables
close all

% Choose basic settings for simulations
n_t = 1e2;
delta_t = 1e-1;
f_fname = @f_lin2D;
g_fname = @g_Id;

u       = zeros(1,n_t);


% Parameters of the simulation
alpha   = 1e1;
sigma   = 1e1;
theta   = 1;
phi     = [];

% Build priors for model inversion
priors.muX0 = zeros(2,1);
priors.SigmaX0 = 1e-0*eye(2);
priors.muTheta = 0*ones(1,1);
priors.SigmaTheta = 1e-1*eye(1);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% render 1st state deterministic
iQx = diag([1e2;1]);
for t=1:n_t
    priors.iQx{t} = iQx;
end


% Build options and dim structures for model inversion
options.priors      = priors;
options.inF.deltat = delta_t;
options.inF.b = 5e-1;
options.backwardLag  = 16;
dim.n_theta         = 1;
dim.n_phi           = 0;
dim.n               = 2;

% options.checkGrads = 1;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,priors.muX0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);

% Call inversion routine
% [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options);
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma);

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(...
        n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end


