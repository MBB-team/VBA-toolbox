function [posterior, out] = demo_stochasticModel ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_stochasticModel ()
% Demo for stochastic system with binomial output
%
% This demo simulates and inverts a model of a stochastic dynamical system, 
% which is observed through a nonlinear sigmoid mapping (binary 
% observations).
%
% /////////////////////////////////////////////////////////////////////////

% number of observations
N = 1e2;

%% Specify the model
% =========================================================================

% AR(1) evolution function
function [fx, dfdx, dfdp] = f_evolution(x, ~, ~, ~)
    fx = x;
    dfdx = eye(length(x));
    dfdp = [];
end

% binary observations from a biased sigmoid mapping
function [gx, dgdx, dgdp] = g_observation(x, P, ~, in)
    [gx, dgdx, dsdp] = VBA_sigmoid(x,'slope',exp (P), in, 'derivatives', {'slope'});
    dgdx = dgdx';
    dgdp = dsdp * exp (P);
end

% bias of the sigmoid mapping
options.inG.center = randn(4,1);

% binary observations
options.sources.type = 1; 

%% Simulate data
% =========================================================================

% precision of stochastic innovation 
alpha   = 1e1;
% parameters
theta = [];
phi = 0.5;
x0 = 0;

% simulate
[y,x,x0,eta,e] = VBA_simulate (N,@f_evolution,@g_observation,theta,phi,[],alpha,[],options,x0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);

%% Inversion
% =========================================================================
% override default priors (deterministic -> stochastic evolution)
options.priors.a_alpha = 1;
options.priors.b_alpha = 1;

% window of the hidden state estimation filtering
options.backwardLag = 5;

% dimensions of the problem
dim.n_theta         = 0;
dim.n_phi           = 1;
dim.n               = 1;

% inversion routine
[posterior, out] = VBA_NLStateSpaceModel(y,[],@f_evolution,@g_observation,dim,options);

%% Results
% =========================================================================
% Display results
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,[]);

end