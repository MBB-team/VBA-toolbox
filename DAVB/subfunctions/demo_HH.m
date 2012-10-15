% Demo for Hodgkin-Huxley PA propagation

clear variables
close all

% Choose basic settings for simulations
n_t = 8e2;
delta_t = 1e-1;         % 1Hz sampling rate
f_fname = @f_HH;
g_fname = @g_Id;

u       = 1e0*randn(1,n_t)>2;
% figure,plot(u)

% Build options structure for temporal integration of SDE
inF.delta_t = delta_t;
inG.ind         = 1;
options.inF     = inF;
options.inG     = inG;
options.decim = 1;


% Parameters of the simulation
alpha   = Inf;
sigma   = 1e4;
phi     = [];
x0 = [0;-3;-3/4;1/3];
theta = [1;0*ones(4,1)];

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
disp('--paused--')
pause

dim.n = 4;
dim.n_phi = 0;
dim.n_theta = 4;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);


% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)


