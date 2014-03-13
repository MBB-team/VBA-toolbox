% this script reproduces a MATLAB statistical toolbox demo

clear all
close all


x = [5.72 4.22 5.72 3.59 5.04 2.66 5.02 3.11 0.13 2.26 ...
     5.39 2.57 1.20 1.82 3.23 5.46 3.15 1.84 0.21 4.29 ...
     4.61 0.36 3.76 1.59 1.87 3.14 2.45 5.36 3.44 3.41]';
y = [2.66 2.91 0.94 4.28 1.76 4.08 1.11 4.33 8.94 5.25 ...
     0.02 3.88 6.43 4.08 4.90 1.33 3.63 5.49 7.23 0.88 ...
     3.08 8.12 1.22 4.24 6.21 5.48 4.89 2.30 4.13 2.17]';

priors.muPhi = [0;0];
priors.SigmaPhi = 1e2.*eye(2);
priors.a_sigma = 1;
priors.b_sigma = 1;
options.priors = priors;
options.inG.x = x;
dim.n = 0;
dim.n_phi = 2;
dim.n_theta = 0;
g_fname = @g_nl0;
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);