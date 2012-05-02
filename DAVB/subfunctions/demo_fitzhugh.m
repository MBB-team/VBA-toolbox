% FitzHugh-Nagum demo for calcium imaging of spike trains

clear variables
close all

% Choose basic settings for simulations
n_t = 4e2;
deltat = 1e-1;         % 10Hz sampling rate
f_fname = @f_FitzHughNagumo;
g_fname = @g_Id;

u       = 0e0*(randn(1,n_t)>2);
u = randn(1,n_t);
figure,plot(u)

% Build options structure for temporal integration of SDE
inF.delta_t     = deltat;
inF.a           = 0.5;
inG.ind         = 1;
options.inF     = inF;
options.inG     = inG;
options.decim   = 10;


% Parameters of the simulation
alpha   = Inf;
sigma   = Inf;
theta   = [1;0;0*randn(3,1)];
phi     = [];



% Build priors for model inversion
priors.muX0 = 0*ones(3,1);
priors.SigmaX0 = 1e-1*eye(3);
priors.muTheta = theta;%0.*ones(5,1);
priors.SigmaTheta = 1e1*eye(5);
priors.SigmaTheta(2,2) = 0;
priors.a_alpha = Inf;
priors.b_alpha = 0;
priors.a_sigma = 1e1;
priors.b_sigma = 1e-4;

options.priors      = priors;
dim.n_theta         = 5;
dim.n_phi           = 0;
dim.n               = 3;

% Build time series of hidden states and observations
x0 = [0;0.8;-0.7];
[y,x,x0,eta,e] = simulateNLSS(...
    n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
disp('--paused--')
pause

priors.a_alpha      = 1e0;
priors.b_alpha      = 1e0;
priors.a_sigma      = 1e3;
priors.b_sigma      = 1e0;
for t = 1:n_t
    dq              = 1e4*ones(3,1);
    dq(2)           = 1;
    priors.iQx{t}   = diag(dq);
end
options.priors = priors;


% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],f_fname,g_fname,dim,options);


%------------ Display results ------------------%
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)

