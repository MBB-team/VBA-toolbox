% Demo for linear system with AR(1) state noise
% The class of generative models that the toolbox can deal with is actually
% more general that may appear at first glance. In fact, it can deal with
% any form of of auto-regressive or state-dependant noises (both at the
% level of hidden states and observations). In most instances, it suffices
% to construct an augmented state space X, where Xt = (xt,zt)T and
% appropriately modify the evolution and observation functions, as well as
% the priors. In the case of AR(1) noise, this could be implemented as
% follows: 
%      X_t+1 = [ f(x_t,theta,u_t) + zt , zt ] + eta_t
%      yt    = g(x_t,phi) + e_t
% where the f is the original evolution function. Note that one has to
% ensure that z is the most likely perturbation force on the system. This
% is because both AR(1) and white noise (respectively zt and ?t) can drive 
% the system. This can be done by defining the augmented state noise
% covariance matrix Qx (of eta) such that its upper-left half is
% effectively zero. In addition, one may have to increase the lag k. This
% is because the effect of the AR(1) retarded state noise on the hidden
% states is maximal one time step in the future. On thus need to look back
% one time step in the past to infer on the retarded state noise.

clear variables
close all

% Choose basic settings for simulations
n_t = 2e2;
dt = 1e-1;
f_fname = @f_embedAR;  % this is an AR(1) embedding evolution function
g_fname = @g_embedAR;  % this is an AR(1) embedding observation function
u = zeros(1,n_t);
in.opt.f_fname = @f_lin2D; % this is the native evolution function (n=2,n_theta=1)
in.opt.g_fname = @g_Id;
in.opt.priors.a_alpha = 1;
in.opt.priors.b_alpha = 1;
in.dim.n = 2;
in.dim.n_theta = 1;
in.dim.n_phi = 0;
in.dim.p = 2;
in.dim.n_t = n_t;
[in.opt] = VBA_check([],u,in.opt.f_fname,in.opt.g_fname,in.dim,in.opt);
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
options.backwardLag  = 8;
dim.n_theta = in.dim.n_theta;
dim.n_phi = 0;
dim.n = 2*in.dim.n;
dim.p = 2;
dim.n_t = n_t;
options.dim = dim;


% Build time series of hidden states and observations
% [y,x,x0,eta,e] = simulateNLSS_old(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)


% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma)


% now invert data without AR(1) priors on state noise
f_fname = @f_lin2D;
g_fname = @g_Id;
dim.n = in.dim.n;
options.inF = in.opt.inF;
options.inG = [];
options.priors.muX0 = zeros(in.dim.n,1);
options.priors.SigmaX0 = 1e-1*eye(in.dim.n);
for t=1:n_t
    options.priors.iQx{t} = eye(2);
end
[p0,o0] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(p0,o0,y-e,x(1:2,:),x0(1:2),theta,phi,alpha,sigma)




