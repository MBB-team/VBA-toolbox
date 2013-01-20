function [] = demo_calciumImaging()

% Demo for calcium imaging of spike trains

clear variables
close all

% Choose basic settings for simulations
% n_t = 8e2;
n_t = 1e3;

delta_t = 2e-3;         % 10Hz sampling rate
% delta_t = 1e0;         % 10Hz sampling rate

% f_fname = @f_lin1D;
f_fname = @f_HH_calcium;
g_fname = @g_Id;

% u       = 1e0*randn(1,n_t)>2;
u       = 5e1*(randn(1,n_t)>1.5);

figure,plot(u)

% Build options structure for temporal integration of SDE
inF.delta_t = delta_t;
inF.sc = 5e1; % time scale change for HH model
inG.ind         = 1;
options.inF     = inF;
options.inG     = inG;
options.decim = 1;


% Parameters of the simulation
alpha   = Inf;
sigma   = 1e2;
phi     = [];
x0 = [0;0;-3;-3/4;1/3];
theta = [1;1;0*ones(3,1)];

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
% 
% dim.n = 5;
% dim.n_theta = 5;
% dim.n_phi = 0;
% options.priors.a_alpha = Inf;
% options.priors.b_alpha = 0;
% options.priors.muX0 = x0;
% [posterior,out] = VBA_NLStateSpaceModel(y(:,1:500),u(:,1:500),f_fname,g_fname,dim,options);
% 
return
% 
% disp('--paused--')
% pause

Y.calcium = y;
Y.calcium_time = [0;inF.delta_t];
[posterior,out] = detectSpikes(Y,'HH',12);



