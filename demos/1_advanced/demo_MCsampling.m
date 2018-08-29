% this script demonstrates how to derive and display Monte-Carlo esimates
% of the (prior) predictive density over hidden states of dynamical
% systems.

clear variables
close all

% Choose basic settings for simulations
n_t = 5e2;
deltat = 2e-2;
alpha   = 1e0/deltat;
sigma   = 1e-1;
theta   = [28;10;8/3];
f_fname = @f_Lorenz;
g_fname = @g_sigmoid;
u       = [];


% Build options structure for temporal integration of SDE
inG.scale = 50;
inG.slope = 0.2;

inF.deltat = deltat;
options.inF     = inF;
options.inG     = inG;

% Build priors
priors.muX0 = 1*ones(3,1);
priors.SigmaX0 = 1e0*eye(3);
priors.muTheta = [28;10;8/3];
priors.SigmaTheta = eye(3);
priors.a_alpha = 1/deltat;
priors.b_alpha = 1;
priors.a_sigma = 1;
priors.b_sigma = 1;

% Build options and dim structures for model inversion
options.priors = priors;
dim.n_theta = 3;
dim.n_phi  = 0;
dim.n = 3;
dim.p = dim.n;

[pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim,1e2);

[h] = plotDensity(f_fname,g_fname,u,n_t,options,dim,pX,gX,pY,gY);

