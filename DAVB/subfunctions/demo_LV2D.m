% Demo for Lotak-Volterra competitive system

clear variables
close all

% Choose basic settings for simulations
f_fname = @f_LV2D;
g_fname = @g_Id;
u       = [];
n_t = 1e2;
deltat = 1e-1;
alpha   = 5e2;
sigma   = 1e1;
theta   = 2e-1*[1;1;1;1];
phi     = [];


% Build options structure
priors.muX0 = [theta(1)./theta(2);theta(3)./theta(4)+1];
priors.SigmaX0 = 1e-4*eye(2);
priors.muTheta = zeros(length(theta),1);
priors.SigmaTheta = 1e0*eye(4);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;
options.priors = priors;
options.decim = 10;
inF.deltat = deltat;
options.inF     = inF;
options.backwardLag = 1;

% Build options and dim structures for model inversion
dim.n_theta = 4;
dim.n_phi   = 0;
dim.n       = 2;



% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
% disp('--paused--')
% pause

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




