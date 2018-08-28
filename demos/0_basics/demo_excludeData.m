function demo_excludeData ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_excludeData ()
% Demo of simulation and inversion of a dynamical system that is sampled on
% an irregular grid (not all timepoints are observed)
%
% More generally, this demo shows how to deal with data exclusion
%
% /////////////////////////////////////////////////////////////////////////

% number of simulated points
N = 257;

%% Specify the model
% =========================================================================
% For the sake of the demonstration, we use a simple linear dynamical
% system $ dx/dt = Ax + phi' u $ that is directly observed $ y = x + noise $

% evolution function
% -------------------------------------------------------------------------
options.inF.A = [- 4, - 16; 4, - 4];
options.inF.dt = 1e-1;

function fx = f_evolution(x,P,u,in)
    xdot = in.A * x + diag(P) * u;
    fx = x + in.dt * xdot;    
end

f_fname = @f_evolution;

% observation function
% -------------------------------------------------------------------------
g_fname = @g_Id;

%% Simulate data
% =========================================================================

% inputs
u = randn(2, N);

% parameters
x0 = [0; 0]; % initial state
theta = [1; 2]; % effect of inputs
phi = []; % no observation parameters
alpha = Inf; % deterministic
sigma = 1; % observation noise

% Build full time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (N,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

% display full time series of hidden states and observations
displaySimulations(y,x,eta,e);

%% Decimation
% =========================================================================
% here we will simulate a sampling of data on an irregular grid  by 
% 'removing' data from the full simulation

% sampling grid (exponential timestps)
timepoints = 2 .^ (1 : floor (log (N) ./ log (2)));

% exclude all data but timepoints: 
% if options.isYout = 1, the datapoint is ignored
options.isYout = ones (size (y));
options.isYout(:, timepoints) = 0;

%% Inversion
% =========================================================================

% dimensions of the problem
dim.n_theta = 2;
dim.n_phi = 0;
dim.n = 2;
dim.p = 2;

% Invert deterministic model
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% display
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma);

end