% this demo replicates the exemple in the 'dynamical system theory' part
% of the 'principles of DCM' talk (SPM course, may 2013).

% clear all
close all

nt = 1e3;
f_fname = @f_2d;
g_fname = @g_Id;

in.A = [-2 -16
        4   -2];
in.dt = 1e-2;
x0 = [1;1];
theta = [];
phi = [];
alpha = Inf;
sigma = Inf;
u =zeros(1,nt);

% Build options and dim structures for model inversion
options.inF = in;
options.inG = in;
dim.n_theta         = 0;
dim.n_phi           = 0;
dim.n               = 2;


% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(nt,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);



% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
