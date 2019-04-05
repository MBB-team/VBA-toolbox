% demo for calcium imaging of spike trains
% The script first simulates the response of a Fitz-Hugh-Nagumo (FHN)
% neuron to spiky input current. It then inverts a FHN-neuron model,
% without the input current info. Practically speaking, this means
% deconvolving the FHN-neuron response to estimate its input.
% [see demo_HodgkinHuxley.m and demo_FHN.m]

clear variables
close all

% Choose basic settings for simulations
n_t = 1e2;
deltat = 1e-1;         % 10Hz sampling rate
f_fname = @f_FitzHughNagumo;
g_fname = @g_Id;

u = randn(1,n_t);

% Build options structure for temporal integration of SDE
inF.dt     = deltat;
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

% Build time series of hidden states and observations
x0 = [0;0.8];
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause

% Build priors for model inversion
priors.muX0 = 0*ones(2,1);
priors.SigmaX0 = 1e0*eye(2);
priors.muTheta = 0.*ones(5,1);
priors.SigmaTheta = 1e0*eye(5);
priors.SigmaTheta(2,2) = 0;
priors.a_alpha      = 1e0;
priors.b_alpha      = 1e0;
priors.a_sigma      = 1e0;
priors.b_sigma      = 1e0;
for t = 1:n_t
    dq              = 1e4*ones(2,1);
    dq(1)           = 1;
    priors.iQx{t}   = diag(dq);
end
options.priors = priors;

dim.n_theta         = 5;
dim.n_phi           = 0;
dim.n               = 2;


% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],f_fname,g_fname,dim,options);


%------------ Display results ------------------%
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma);


[ehat,v_e,etahat,v_eta] = VBA_getNoise(posterior,out);
[haf,hf,hp] = plotUncertainTimeSeries(etahat,VBA_getVar(v_eta));
set(get(haf,'parent'),'name','simulated versus estimated input')
set(haf,'nextplot','add')
rescale = inv(u*u')*u*etahat(1,:)';
plot(rescale*u,'k--')



