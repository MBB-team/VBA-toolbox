function [posterior, out] = demo_noiseAR1 ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_noiseAR1 ()
% Demo for linear system with auto-regressive AR(1) state noise
%
% The class of generative models that the toolbox can deal with is actually
% more general that may appear at first glance. In fact, it can deal with
% any form of of auto-regressive or state-dependant noises (both at the
% level of hidden states and observations). In most instances, it suffices
% to construct an augmented state space X, where Xt = (xt,zt) and
% appropriately modify the evolution and observation functions, as well as
% the priors. In the case of AR(1) noise, this could be implemented as
% follows: 
% $$
%      X_t+1 = [ f(x_t,theta,u_t) + z_t ; z_t ] + eta_t
%      y_t     = g(x_t,phi) + e_t
% $$
% where the f is the original evolution function.
% Note that both AR(1) and white noise (respectively z_t and eta_t) can 
% drive the system. To ensure that z is the most likely perturbation force
% on the system, one can define the augmented state noise covariance matrix
% Qx (of eta) such that its upper-left half is effectively zero.
% In addition, one may have to increase the lag k of the  estimation 
% procedure. This is because the effect of the AR(1) retarded state noise 
% on the hidden states is maximal one time step in the future. One thus 
% need to look back one time step in the past to infer on the retarded
% state noise.
%
% /////////////////////////////////////////////////////////////////////////

% number of observations
N = 1e2;

% timestep
dt = 1e-1;

%% Specify the model
% =========================================================================
% we use the generic embeding functions for AR(1) noise: the original 
% evolution and observation function are passed as options to the embeding
% scheme model

% embedding function
% -------------------------------------------------------------------------
f_fname = @f_embedAR;  % this is an AR(1) embedding evolution function
g_fname = @g_embedAR;  % this is an AR(1) embedding observation function

% definition of the native model
% -------------------------------------------------------------------------

% evolution and observation functions
in.opt.f_fname = @f_lin2D; 
in.opt.g_fname = @g_Id;

% options structures
in.opt.inF.deltat = dt;
in.opt.inF.b = 5e-1;

% dimension
in.dim.n = 2;
in.dim.n_theta = 1;
in.dim.n_phi = 0;
in.dim.p = 2;
in.dim.n_t = N;

% trigger the initialization of priors and options for a stochastic system
u = zeros(1,N);
in.opt.priors.a_alpha = 1;
in.opt.priors.b_alpha = 1;
in.opt = VBA_check([],u,in.opt.f_fname,in.opt.g_fname,in.dim,in.opt);

%% Simulate data
% =========================================================================

% Parameters of the simulation
x0 = cat(1,[0; 0], [0; 0]); % [x_t; z_t]
theta = 1;
phi = [];
alpha = 1e1;
sigma = 1e-1;

% pass on option sructures
options.inF = in;
options.inG = in;

% % increase precision on non-AR(1) state noise component and therefore 
% % increase the weight of z_t on system innovation
% % -------------------------------------------------------------------------
% iQxt = 1e2 * eye (in.dim.n); % precision of eta_t component added to x_t
% iQzt = eye (in.dim.n); % precision of eta_t component added to z_t
% iQx = blkdiag (iQxt, iQzt);
% 
% for t = 1 : N
%     options.priors.iQx{t} = iQx;
% end

% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (N,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);

%% Inversion: accounting for AR(1) noise
% =========================================================================

% dimension of the problem
dim = in.dim;
dim.n = 2 * dim.n; % [x_t; z_t]

% stochastic evolution
options.priors.a_alpha = 1;
options.priors.b_alpha = 1;

% shrinkage prior on X0
options.priors.SigmaX0 = 1e-1 * eye(dim.n);

% lag of the inversion scheme
options.backwardLag  = 16;

% Call inversion routine with AR(1) priors on state noise
[posterior.AR, out.AR] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
[h(1), h(2)] = displayResults(posterior.AR,out.AR,y-e,x,x0,theta,phi,alpha,sigma);
h(3) = out.AR.options.hf;
set (h, 'Name', 'Results: AR(1) noise');


%% Inversion: white noise only
% =========================================================================
% now we invert data without AR(1) priors on state noise

% use native dynamics
f_fname = @f_lin2D;
g_fname = @g_Id;
% no state expansion
dim.n = in.dim.n;
% undo embedding
options.inF = in.opt.inF;
options.inG = in.opt.inG;
% reset priors on state innovation
for t = 1 : N
    options.priors.iQx{t} = eye(2);
end
% shrinkage prior on X0
options.priors.SigmaX0 = 1e-1 * eye(dim.n);

% Call classical inversion routine (white noise)
[posterior.WN,out.WN] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
[h(1), h(2)] = displayResults(posterior.WN,out.WN,y-e,x(1:2,:),x0(1:2),theta,phi,alpha,sigma);
h(3) = out.WN.options.hf;
set (h, 'Name', 'Results: white noise');

