% Demo for partially observable oscillatory system.
% This demo inverts a model of a linear oscillatory system, which is
% observed through a nonlinear sigmoid function.

clear variables
close all

% Choose basic settings for simulations
n_t = 2e2;
delta_t = 2e-1;         % integration time step (Euler method)
f_fname = @f_lin2D;
g_fname = @g_sigmoid;

u       = zeros(1,n_t);

% Build options structure for temporal integration of SDE
inF.deltat      = delta_t;
inF.a           = 0.1;
inF.b           = 0.9e-2;
inG.scale = 2;
inG.y0          = -1;
inG.ind         = 1; % only x(1) is partially observable 
options.inF     = inF;
options.inG     = inG;


% Parameters of the simulation
alpha   = 1e1;
sigma   = 1e1;
theta   = 1;
phi     = [];

% Build priors for model inversion
priors.muX0 = zeros(2,1);
priors.SigmaX0 = 1e0*eye(2);
priors.muTheta = 0*ones(1,1);
priors.SigmaTheta = 1e0*eye(1);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

for t=1:n_t
    priors.iQx{t} = eye(2);
    priors.iQx{t}(1,1) = 1e2;
end

% Build options and dim structures for model inversion
options.priors      = priors;
options.backwardLag  = 16;
dim.n_theta         = 1;
dim.n_phi           = 0;
dim.n               = 2;

% options.checkGrads = 1;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause

% Call inversion routine
% [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options);
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(...
        n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end


