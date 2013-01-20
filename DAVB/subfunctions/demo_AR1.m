% Demo for linear system with AR(1) state noise

clear variables
close all

% Choose basic settings for simulations
n_t = 2e2;
dt = 1e-1;
f_fname = @f_embedAR;  % this is an AR(1) embedding evolution function
g_fname = @g_embedAR;  % this is an AR(1) embedding observation function
u = [];
in.opt.f_fname = @f_lin2D; % this is the native evolution function (n=2,n_theta=1)
in.opt.g_fname = @g_Id;
in.dim.n = 2;
in.dim.n_theta = 1;
in.dim.n_phi = 0;
in.dim.p = 2;
in.dim.n_t = n_t;
[in.opt] = VBA_check([],u,in.opt.f_fname,in.opt.g_fname,in.dim,[]);
in.opt.inF.deltat = dt;
in.opt.inF.b = 5e-1;

% Parameters of the simulation
alpha   = 1e1;
sigma   = 1e-1;
theta   = 1;
phi     = [];

% Build priors for model inversion
priors.muX0 = zeros(2*in.dim.n,1);
priors.SigmaX0 = 1e-1*eye(2*in.dim.n);
priors.muTheta = zeros(in.dim.n_theta,1);
priors.SigmaTheta = 1e0*eye(in.dim.n_theta);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;
% increase precision on non-AR(1) state noise component
iQx = ones(2*in.dim.n,1);
iQx(1:in.dim.n) = 1e2;
iQx = diag(iQx);
for t=1:n_t
    priors.iQx{t} = iQx;
end

% Build options and dim structures for model inversion
options.priors = priors;
options.inF = in;
options.inG = in;
options.backwardLag  = 0;
dim.n_theta = in.dim.n_theta;
dim.n_phi = 0;
dim.n = 2*in.dim.n;
dim.p = 2;
dim.n_t = n_t;
options.dim = dim;


% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)


% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma)



