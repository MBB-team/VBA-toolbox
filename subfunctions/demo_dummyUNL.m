% demo: VBA inversion of generative models with unnormalized likelihoods.
% This demo considers a dummy scenario, whereby a Gaussian log-likelihood
% is specified as if it was arbitrary, i.e. in terms of a quadratic error
% function. We can then compare the VBA_UNL inversion of this model and its
% normalized variant, directly specified in terms of the observation
% function.

close all
% clear all

n = 32;
phi = [1;1;1];
sigma = 1e4;
e = randn(1,n);
u = randn(2,n);
u(2,:) = u(1,:) + u(2,:); % induce parameter non-identifiability
y = phi(1)+ u(1,:)*phi(2)+ u(2,:)*phi(3) + (1/sqrt(sigma))*e;

g_fname = @U_dummy;
dim.n_phi = 3;
options = [];
options.priors.a_sigma = 1/var(y);
options.priors.b_sigma = 1e0;
options.priors.SigmaPhi = 1e0*eye(3);

options.checkGrads = 0;
tic,[p1,o1] = VBA_UNL0(y,u,g_fname,dim,options);toc
set(o1.options.hf,'tag','0');
displayResults(p1,o1,y,[],[],[],phi,[],sigma)

dim.n = 0;
dim.n_theta = 0;
g_fname = @g_Udummy;
[p2,o2] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);
set(o2.options.hf,'tag','0');
