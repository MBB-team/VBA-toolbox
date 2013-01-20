% Demo for Double Well system.
% This demo inverts a model of chaotic double-Well system, which is
% observed through a nonlinear sigmoid observation function.

clear variables
close all

% Choose basic settings for simulations
f_fname = @f_doubleWell;
g_fname = @g_sigmoid;
u       = [];
n_t = 2e2;
deltat = 2e-2;
alpha   = 6e-1/deltat;
sigma   = 1e2;
theta   = [-2;3;1.5];
phi     = [1;1];


% Build options structure for temporal integration of SDE
inG.G0 = 50;
inG.beta = 0.5;
inF.deltat = deltat;
inF.a   = -2;
inF.b   = 3;
options.inF     = inF;
options.inG     = inG;
options.backwardLag = 5;


% Build priors for model inversion
priors.muX0 = [0;0];
priors.SigmaX0 = 1e-3*eye(2);
priors.muTheta = 0.*ones(length(theta),1);
priors.muTheta(3) = 1.5;
priors.SigmaTheta = 1e-1*eye(3);
% priors.SigmaTheta(3,3) = 0;
priors.muPhi = 1*ones(length(phi),1);
priors.SigmaPhi = 0e0*eye(2);
priors.SigmaPhi(1,1) = 0;
priors.a_alpha = 1e2;
priors.b_alpha = 1e2;
priors.a_sigma = 1e2;
priors.b_sigma = 1e2;

% priors.iQx = cell(n_t,1);
% for t=1:n_t
%     priors.iQx{t} = eye(2);
%     priors.iQx{t}(2,2) = 1e1;
% end

% Build options and dim structures for model inversion
options.priors = priors;
dim.n_theta = 3;
dim.n_phi   = 2;
dim.n       = 2;



% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
% disp('--paused--')
% pause

% Call inversion routine
% options.EvoEval = 'show_potential(posterior)';
% options.DisplayWin = 0;
% options.gradF = 1;
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




