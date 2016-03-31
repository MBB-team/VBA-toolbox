close all
clear all

dim.p = 32;
dim.n_t = 1;
dim.n_phi = 16;
dim.n_theta = 0;
dim.n = 0;


g_fname = @g_classif;
options.binomial = 1;
options.priors.muPhi = zeros(dim.n_phi,1);
options.priors.SigmaPhi = 1e0*eye(dim.n_phi);
options.isYout = zeros(dim.p,1);
options.DisplayWin = 1;
options.verbose = 0;

options.inG.X = randn(dim.n_phi-1,dim.p);
phi = randn(dim.n_phi,1);
[y,x,x0,eta,e] = simulateNLSS(dim.n_t,[],g_fname,[],phi,[],[],[],options,[]);
g = y-e;
g = g>0.5; % denoised data

[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);