% Demo for Fitz-Hugh-Nagumo PA propagation
% Short bursts of depolarizing input current (u) are sent to a neuron that
% responds according to Fitz-Hugh-Nagumo model. Note that this does not
% generates 'all-or-none' spikes (as in the Hodgkin-Huxley model), but
% supra-threshold input create large deviations from equilibrium membrane
% depolarization. NB: an AP approximately corresponds to a deviation of
% about 1 A.U. w.r.t. equilibrium.
% [see demo_fitzhugh.m]

clear variables
close all



% Choose basic settings for simulations
n_t = 1e3; % # length of time series 
dt = 1e-1; % 10Hz sampling rate
f_fname = @f_FitzHughNagumo;
g_fname = @g_Id;

u = 1e1*(randn(1,n_t)>2); % about 10 to 15% of samples are 'on'

% Build options structure for temporal integration of SDE
inF.dt = dt;
inG.ind         = 1;
options.inF     = inF;
options.inG     = inG;
options.decim = 1;


% Parameters of the simulation
alpha   = Inf;
sigma   = 1e2;
phi     = [];
x0 = [1.2;0.62];
theta = [0.4;0.1*ones(3,1)];

% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);


hf = figure('color',[1 1 1]);
ha = subplot(3,1,1,'parent',hf,'nextplot','add');
plot(ha,u')
title(ha,'input current')
xlabel('time')
ha = subplot(3,1,2,'parent',hf,'nextplot','add');
title(ha,'recovery variable')
plot(ha,x(2,:))
xlabel('time')
ha = subplot(3,1,3,'parent',hf,'nextplot','add');
plot(ha,x(1,:))
title(ha,'output membrane depolarization (A.U.)')
xlabel('time')


% invert FHN model on membrane potential data
dim.n = 2;
dim.n_phi = 0;
dim.n_theta = 4;
priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);


% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma);


