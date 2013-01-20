%% Demo_QLearning_2Q.m
%%
% This script operates the inversion of a reinforcement learning model with
% a softmax decision rule.
% TASK DESCRIPTION:
% At each trial, the agent chooses between two alternative actions, which
% are rewarded according to a predefined (probabilistic) reward schedule.
% This is a typical two-armed bandit task.
% MODEL DESCRIPTION:
% The model is a typical Rescorla Wagner learning rule, whereby the action
% value is updated proportionaly to the reward prediction error, i.e.:
% Q(t+1) = Q(t) + alpha*(R(t)-Q(t))
% where Q(t) is the learned value at trial t, alpha is the learning rate,
% and R(t) is the delivered reward. Note that there are two Q-values (one
% for each action). The decision is emitted according to a softmax rule.

clear all 
close all
clc

%% Simulate data

% 1- Q-learner parameters
Ntrials = 200; % # trials
alpha = 0.2; % learning rate
beta = 1; % inverse temperature
x0 = 1e-1*randn(2,1); % initial Q-values

% 2- Define feedback.
% NB: Here, the feedback is defined through a feedback function, which
% delivers a reward R1 (resp. R2) to the agent when she chooses action 1
% (resp. action 2). This is a (very simple) special case of a more general
% procedure, by which feedback could be delivered according to a complex
% rule.
% NB2: this is a probabilistic reward schedule, and the reward
% contingencies change at half of the session. 
R = rand(2,Ntrials);
R(1,1:Ntrials/2) = R(1,1:Ntrials/2)<0.8;
R(1,Ntrials/2+1:end) = R(1,Ntrials/2+1:end)<=0.2;
R(2,1:Ntrials/2) = R(2,1:Ntrials/2)<=0.2;
R(2,Ntrials/2+1:end) = R(2,Ntrials/2+1:end)<0.8;

% 3- options for data simulations
fb.inH.u0 = R; % reward scehdule
fb.h_fname = @h_reward_2Q; % feedback function
fb.indy = 1; % index of chosen action
fb.indfb = 2; % index of received feedback
f_fname = @f_Qlearn_2Q;
g_fname = @g_softmax_2Q;
theta = sigm(alpha,struct('INV',1));
phi = [log(beta);0.5];
u = [zeros(max([fb.indfb(:);fb.indy(:)]),Ntrials)];
alpha = Inf;
sigma =  Inf;
in.param_transform.type = 'modified sigmoid';
in.param_transform.a = 0;
in.param_transform.b = 1;
options.inF = in;
options.inG = in;
options.binomial = 1; % Dealing with binary data
options.isYout = zeros(1,Ntrials); % Excluding data points
dim = struct('n',2,...  % hidden states = Q-values
    'p',1,... % output = emitted choices
    'n_theta',1,... % evolution parameters
    'n_phi', 2,... % observation parameters
    'n_t',Ntrials);
options.dim = dim;
options.skipf = zeros(1,Ntrials);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

% 4- simulate and display data
[y,x,x0,eta,e,u] = simulateNLSS_fb(Ntrials,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);
displaySimulations(y,x,eta,e)


%% Invert model

% 1- Define priors (i.i.d. normal densities)
priors.muPhi = zeros(dim.n_phi,1); 
priors.SigmaPhi = 1e1*eye(dim.n_phi);
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaTheta = 1e1*eye(dim.n_theta);
priors.muX0 = ones(dim.n,1)*0;
priors.SigmaX0 =1e1*eye(dim.n);
% priors.a_alpha = Inf; % No state noise for deterministic update rules
% priors.b_alpha = 0;
options.priors = priors;

% 2- call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%% evaluate inversion results
displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)
