% This script operates the inversion of a temporal discounting model with
% a softmax decision rule.
% TASK DESCRIPTION:
% At each trial, two alternatives are proposed to the agent, each of which
% has two dimensions (reward value and delay).
% MODEL DESCRIPTION:
% Agents are assumed to devaluate items accroding to a hyperbolic discount
% rate, i.e.:
% V(d) = V0 / (1 + k*d),
% where d is the delay, V0 is the reward magnitude, and k is the delay
% discounting parameter. The decision is emitted according to a softmax
% rule.

close all
clear variables
clc

% 1- agent parameters
K = 1/15; % discount parameter
beta = 0.5; % inverse temperature
N = 100; % # trials

% 3- options for data simulations
% 3.1 evolution/observation functions
g_fname = @g_1Dhyp; % observation function, Parameters : [K, beta] => n_phi = 2
f_fname = []; % no learning
% 3.2 define set of alternative choices
T = zeros(2,N); % reward delays
T(2,:) = rand(1,N)*30;
V = zeros(2,N); % reward values
V(1,:) = 10;
V(2,:) = 20;
inG.T = T; % 2 * N (times of reception)
inG.V = V;  % 2 * N  (objective values)
% 3.3 OPTIONS structure
dim = struct('n',0,'n_theta',0,'n_phi',2,'p',N,'n_t',1);
options.dim = dim;
options.inG = inG; 
options.binomial = 1; % Dealing with binary data
options.verbose = 0;
options.inG = inG; 
% 3.4 simulate data
[y,x,x0,eta,e] = simulateNLSS(N,f_fname,g_fname,[],[log(K),log(beta)],[],Inf,Inf,options,[]);

% 4- Invert model
% 4.1 Define priors (i.i.d. normal densities)
priors.muPhi = zeros(dim.n_phi,1); 
priors.SigmaPhi = 1e0*eye(dim.n_phi);
options.priors = priors;
% 4.2 call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],f_fname,g_fname,dim,options);
% 4.3 evaluate inversion results
displayResults(posterior,out,y,[],[],[],[log(K),log(beta)],Inf,Inf)

    
    
