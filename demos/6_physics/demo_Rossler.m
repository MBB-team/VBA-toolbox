% Demo for Rossler oscillator.
% This demo inverts a model of chaotic Rossler system, which is
% observed through a nonlinear sigmoid observation function.

clear variables
% close all

% Choose basic settings for simulations
f_fname = @f_Rossler;
g_fname = @g_sigmoid;
u       = [];
n_t     = 500;
deltat  = 5e-2;
alpha   = 1e3;
sigma   = 1e3;
theta   = [0.2;0.2;2.791];
phi     = [];


% Build options structure for temporal integration of SDE
inF.deltat = deltat;
options.inF     = inF;
options.inG     = struct;
options.backwardLag = 4;

% Build priors for model inversion
priors.muX0 = 0*[1;1;1];
priors.SigmaX0 = 1e0*eye(3);
priors.muTheta = 0.*ones(length(theta),1);
priors.SigmaTheta = 1e0*eye(3);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% Build options and dim structures for model inversion
options.priors = priors;
dim.n_theta = 3;
dim.n_phi   = 0;
dim.n       = 3;

% Build time series of hidden states and observations
stop = 0;
it = 1;
itmax = 15;
while ~stop
    [y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
    if ~ VBA_isWeird ({x, y}) || it >= itmax
        stop = 1;
    else
        it = it+1;
    end
end

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
VBA_getSubplots ();
% disp('--paused--')
% pause


% Call inversion routine
if ~ VBA_isWeird ({x, y}) 
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
end




