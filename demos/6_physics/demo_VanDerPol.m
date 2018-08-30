% Demo for DAVB Van Der Pol oscillator.
% This demo inverts a model of nonlinear Van Dr Pol oscillator, which is
% observed through a nonlinear sigmoid function.

clear variables
close all

% Choose basic settings for simulations
n_t = 2e2;
deltat = 1e-1;
f_fname = @f_vanDerPol;
g_fname = @g_sigmoid;
alpha   = 1e2;
sigma   = 1e1;
theta   = [1e0];
phi     = [];
u       = [];

% Build options structure for temporal integration of SDE
inG.scale = 50;
inG.slope = 5;
inF.deltat = deltat;
options.inF     = inF;
options.inG     = inG;
options.backwardLag = 5;


% Build priors for model inversion
priors.muX0 = 1e-1*ones(2,1);
priors.SigmaX0 = 1e0*eye(2);
priors.muTheta = 0.*ones(1,1);
priors.SigmaTheta = 1e0*eye(1);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% Build options and dim structures for model inversion
options.priors      = priors;
dim.n_theta         = 1;
dim.n_phi           = 0;
dim.n               = 2;



% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);



% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause
    

% [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options);

% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);




% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma);

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(...
        n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end

