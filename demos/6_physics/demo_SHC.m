% Demo for stable heteroclinic channel (SHC).
% This demo inverst a model of hierarchical SHC partly observable. The
% model has no observation/evolution parameters (except initial
% conditions).

clear variables
close all

% Choose basic settings for simulations
n_t = 5e2;
f_fname = @f_SHC;
g_fname = @g_Id;
u       = [];
deltat  = 2e0;
alpha   = 1e3;
sigma   = 1e4;
theta   = [];
phi     = [];

% Build options structure for temporal integration of SDE
inG.G0 = 50;
inG.ind = 1:4;
inF.deltat = deltat;
options.inF     = inF;
options.inG     = inG;


% Build priors for model inversion
priors.muX0 = [
    -5
   -20
   -20
   -10
   -40
   -20
   -10
];
priors.SigmaX0 = 1e0*eye(7);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% Build options and dim stuctures for model inversion
options.priors      = priors;
options.backwardLag = 5;
options.MaxIterInit = 0; % no deterministic update
dim.n_theta         = 0;
dim.n_phi           = 0;
dim.n               = 7;



% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause


% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma);

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(...
        n_t,theta,phi,zeros(size(u)),alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end


