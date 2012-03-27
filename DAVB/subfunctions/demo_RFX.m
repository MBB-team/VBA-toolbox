% demo random effet analysis (2 levels)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

clear variables
close all

% Choose basic settings for simulations
ns = 12; % number of subjects
n = 256; % number of observations per subject 
f_fname = @g_RFX;
g_fname = @g_RFX;
alpha   = 1e0;
sigma   = 1e0;
theta   = [];
phi     = [];
u       = [];

% Build options structure for temporal integration of SDE
inF.ns = ns;
inF.X = [ones(ns,1),zeros(ns,ns-1)]; % 2d-level design matrix
inG.X = randn(n,ns); % 1st-level design matrix
options.inF     = inF;
options.inG     = inG;

% Build priors for model inversion
priors.muX0 = zeros(ns,1);
priors.SigmaX0 = 0*eye(ns);
priors.SigmaX0(1,1) = 1e2; % group mean effect
priors.a_alpha = 1e2;
priors.b_alpha = 1e2;
priors.a_sigma = 1e2;
priors.b_sigma = 1e2;

% Build options and dim stuctures for model inversion
options.priors      = priors;
options.Laplace = 0;
dim.n_theta         = 0;
dim.n_phi           = 0;
dim.n               = ns;



% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(1,f_fname,g_fname,theta,phi,u,alpha,sigma,options);


% Invert model with non zero 2nd-level effect
[p,o] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(p,o,y,x,x0,theta,phi,alpha,sigma)


% Invert model with zero 2nd-level effect (null hypothesis)
options.priors.SigmaX0 = 0.*options.priors.SigmaX0;
[pH0,oH0] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Re-display alternative model inversion results
VBA_ReDisplay(p,o,1);

