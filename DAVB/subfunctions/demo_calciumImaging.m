% Demo for calcium imaging of spike trains
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

clear variables
close all

% Choose basic settings for simulations
n_t = 4e2;
delta_t = 0.7e-1;         % 10Hz sampling rate
f_fname = @f_lin1D;
g_fname = @g_Id;

u       = randn(1,n_t)>2;
% figure,plot(u)

% Build options structure for temporal integration of SDE
inF.delta_t = delta_t;
inF.a           = 0.5;
inG.ind         = 1;
options.inF     = inF;
options.inG     = inG;



% Parameters of the simulation
alpha   = Inf;
sigma   = 1e4;
theta   = [0.5];
phi     = [];

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
disp('--paused--')
pause

Y.calcium = y;
Y.calcium_time = [0;delta_t];
[posterior,out] = detectSpikes(Y,'HH',12);



