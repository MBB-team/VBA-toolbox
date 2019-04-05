function [posterior, out] = demo_VolterraKernels ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% demo_VolterraKernels ()
% demo of Volterra kernel estimation and use
%
% This code demonstrates how to extract Volterra kernels from a
% dynamical system trajectory. In this particular example, 
% the dynamical system is in fact a Rescorla-Wagner learning agent. Here,
% we use Volterra kernel decomposition to evaluate the long term influence
% of an input (feedback) on the immediate and future responses.
%
% See also demo_dynLearningRate
%
% /////////////////////////////////////////////////////////////////////////

%% Global values
% =========================================================================
N = 1e3;

%% Definition of the model
% =========================================================================
% simple Rescorla-Wagner model
f_fname = @f_rwl;
g_fname = @g_Id;

%% Simulation
% =========================================================================
% parameters
x0 = nan; % random initial state
theta = 0.2; % learning rate
phi = [];
alpha = Inf;
sigma = 1e2;

% options
options = struct ();

% feedback: random positive of negative feedbacks
fb.h_fname = @h_Id; 
fb.inH.u = 2 * VBA_random ('Bernoulli', 0.5, 1, N) - 1;
fb.indy = [];
fb.indfb = 1; % where to store feedbacks in u

% inputs will be fed in by the feedback
u = nan(1, N);

% simulation routine
[y,x,x0,eta,e,u] = VBA_simulate (N,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0,fb);

%% Estimation
% =========================================================================
% dimensions of the problem
dim.n = 1;
dim.n_theta = 1;
dim.n_phi = 0;

% inversion routine
[posterior, out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%% Volterra analysis
% =========================================================================
% compute kernels (NB: you could also find them in out.diagnostics)
maxLag = 16;
kernels = VBA_getVolterraKernels (posterior, out, maxLag);

% compute influence of feedbacks on responses
theoretical = theta * (1 - theta) .^ (0 : maxLag);
empirical.mean = kernels.y.m';
empirical.var = kernels.y.v';

% display
f = VBA_figure('Name', 'Influence of feedback on fugzre responses');
plotUncertainTimeSeries(empirical.mean,sqrt(empirical.var),1, f);
hold on
plot(theoretical,'b')
legend({'estimate (Volterra weight)','credible interval','theoretical'})
box off
xlabel('lag')
ylabel('Influence of input')

