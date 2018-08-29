% this demo replicates the exemple in the 'dynamical system theory' part
% of the 'principles of DCM' talk (SPM course, may 2013).

clear all
close all

nt = 3e2;
f_fname = @f_2d;
g_fname = @g_Id;

in.A = [-2 -16
        4   -2];
in.dt = 1e-2;
x0 = [1;1];
theta = [];
phi = [];
alpha = 1;
sigma = Inf;
u =zeros(1,nt);

% Build options and dim structures for model inversion
options.inF = in;
options.inG = in;
dim.n_theta = 0;
dim.n_phi = 0;
dim.n = 2;
dim.p = 2;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (nt,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);



% display time series of hidden states and observations
displaySimulations(y,x,eta,e);

options.priors.a_alpha = 1;
options.priors.b_alpha = 1;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)