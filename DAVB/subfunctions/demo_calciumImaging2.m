% Demo for calcium imaging of spike trains

clear variables
close all

% Choose basic settings for simulations
% n_t = 8e2;
n_t = 1e2;

delta_t = 1e-1;         % 10Hz sampling rate

f_fname = @f_monoExp;
g_fname = @g_Id;

u       = 1e0*(randn(1,n_t)>1.5);

% figure,plot(u)

% Build options structure for temporal integration of SDE
inF.delta_t = delta_t;
inG.ind         = 1;
options.inF     = inF;
options.inG     = inG;
options.decim = 1;


% Parameters of the simulation
alpha   = 1e5;
sigma   = 1e4;
phi     = [];
x0 = [0];
theta = [1];

% alpha   = Inf;
% sigma   = Inf;
% theta   = [0.5];
% phi     = [];

% x0 = [0];
% theta = [1;0*ones(4,1)];

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)

return
% 
% disp('--paused--')
% pause

Y.calcium = y;
Y.calcium_time = [0;delta_t];
[posterior,out] = detectSpikes(Y,'lognormal',12);



