% Demo for Lorenz choatic system.
% This demo inverts a chaotic Lorenz system model, which is observed
% through a nonlinear sigmoid observation function.

clear variables
close all

% Choose basic settings for simulations
n_t = 3e2;
deltat = 2e-2;
alpha   = 1e0/deltat;
sigma   = 1e-1;
theta   = [28;10;8/3];
f_fname = @f_Lorenz;
g_fname = @g_sigmoid;
phi     = [];
u       = [];


% Build options structure for temporal integration of SDE
inG.scale = 50;
inG.slope = 0.2;

inF.deltat = deltat;
options.inF     = inF;
options.inG     = inG;
% options.u0      = 0*ones(2,1);  % initial condition: input value


% Build priors for model inversion
% for t=1:n_t
%     priors.iQx{t} = eye(3);
%     priors.iQx{t}(1,1) = 1e-3; % makes the 1st state 1000 times more noisy
% end
priors.muX0 = 1*ones(3,1);
priors.SigmaX0 = 1e-1*eye(3);
priors.muTheta = 0.*ones(3,1);
priors.SigmaTheta = 1e1*eye(3);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% Build options and dim structures for model inversion
options.priors      = priors;
options.backwardLag = 1;
dim.n_theta         = 3;
dim.n_phi           = 0;
dim.n               = 3;


options.checkGrads = 0;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause

 
%Call inversion routine
% [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options);
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


