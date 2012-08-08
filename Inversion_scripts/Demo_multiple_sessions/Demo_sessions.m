%-----------------------------------------------------------------------
% This script shows how to deal with multiple independant sessions on an
% operant learning task with a reinforcement learning model (see Demo_Qlearning_2Q.m)
%
% From the generative model of each session is constructed an extended
% model (EM) for all sessions.
%
% Once this model is created, model inversion is similar to that of any
% model
%-----------------------------------------------------------------------

clear all 
close all
clc

Nsessions = 4; % number of sessions
Ntrials = 100; % number of trials

%-----------------------------------------------------------------------
% Preparing feedbacks (task contingencies)
%-----------------------------------------------------------------------
% Probabilistic reward. Correlated probabilities (complementary p1 = 1-p2)
% Fixed probabilities (p1 =0.8) with probability reversal at mid course
Rewards = []; % Matrix of rewards for each alternatives for each session (read it 2 lines by 2)
for i = 1 : Nsessions
R = rand(2,Ntrials);
R(1,1:Ntrials/2) = R(1,1:Ntrials/2)<0.8;R(1,Ntrials/2+1:end) = R(1,Ntrials/2+1:end)<0.2;
R(2,1:Ntrials/2) = R(2,1:Ntrials/2)<0.2;R(2,Ntrials/2+1:end) = R(2,Ntrials/2+1:end)<0.8; % reversal
Rewards = [Rewards;R]; % concatenation through sessions
end

%-----------------------------------------------------------------------
% Model specification (single session)
%-----------------------------------------------------------------------

alpha = 0.1; % learning rate in RL algorithm
beta =3; % inverse temperature in softmax decision rule
U = []; Y = []; IsYout = [];
x0 = [0;0];

f_fname = @f_Qlearn_2Q;
g_fname = @g_softmax_2Q;
theta = sigm(alpha,struct('INV',1));
phi = log(beta);
alpha = Inf; % State noise precision (inf -> deterministic system)

dim = struct('n',2,...  % 2 (Qvalues) * Nsessions
             'p',1,... % total output dimension
             'n_theta',1,... % evolution parameters
             'n_phi', 1,... % observation parameters
             'n_t',Ntrials,... % number of trials
             'u',2); % input data
options.dim = dim; % dimensions of the model

options.inF = []; % input to the evolution function
options.inG = []; % input to the observation function
options.DisplayWin = 1; % display (1) or not (0) the main display window
options.binomial = 1; % Dealing with binary data
options.isYout = zeros(1,Ntrials); % Excluding trials for inversion
                
options.skipf = zeros(1,Ntrials); % Replace evolution function by identity function for trials set to 1;
options.skipf(1) = 1; % apply identity mapping from x0 to x1.


% For any predefined rewards for all alternatives
fb.h_fname = @h_reward_2Q; % handle to put reward of chosen option y into input vector u
fb.indy = 1; % where to put output y in input vector u
fb.indfb = [2]; % indices where put feedbacks in the experimenter data matrix u

%%
%-----------------------------------------------------------------------

%-----------------------------------------------------------------------
%%% SPECIFY MODEL FOR MULTIPLE SESSIONS
%-----------------------------------------------------------------------

% future dimensions of the extended model
dim_e = struct('n_theta',1,... % specify parameter space
             'n_phi', 1);%      

% Information for sessions
in_sessions = struct();
in_sessions.n_sess = Nsessions; % number of sessions
in_sessions.f_fname = @f_Qlearn_2Q; % handle of the shared evolution function
in_sessions.g_fname = @g_softmax_2Q; % handle of the shared observation function
in_sessions.dim_e = dim_e; % specify extended model's parameter space
in_sessions.binomial = 1;
in_sessions.inF = [];
in_sessions.inG = [];
in_sessions.inH = Rewards;
in_sessions.fb = fb;

in_sessions.ind.theta =  ones(Nsessions,1)*[1]; % specify parameter use for each session
in_sessions.ind.phi =   ones(Nsessions,1)*[1]; % specify parameter use for each session

[ f_fname_e,g_fname_e,dim_e,options_e ] = makeExtendedModel(dim,options,in_sessions);



%%
%-----------------------------------------------------------------------
% SIMULATE MODELS
%-----------------------------------------------------------------------
%%%  Option 1 : SIMULATE EACH SESSION INDIVIDUALLY
%-----------------------------------------------------------------------

disp('---------- Simulate each session individually ----------')


for i = 1 : Nsessions
fb.inH.u0 = Rewards((i-1)*2+1:i*2,:); 
u = [zeros(2,Ntrials)];
[y,x,x0,eta,e,u] = simulateNLSS_fb(Ntrials,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);
isYout = zeros(1,Ntrials);
U = [U;u];
Y = [Y;y];
IsYout = [IsYout;isYout];
end

%%
%-----------------------------------------------------------------------
%%%  Option 2 : SIMULATE EXTENDED MODEL
%-----------------------------------------------------------------------

disp('---------- Simulate all sessions together ----------')
x0_e = zeros(dim_e.n,1);
fb_e = options_e.fb;
U_e = zeros(dim_e.u,Ntrials);
for i = 1 : Nsessions
fb_e.inH.sess(i).inH.u0 =  Rewards((i-1)*2+1:i*2,:); % feedback for each session
end

[y,x,x0,eta,e,u] = simulateNLSS_fb(Ntrials,f_fname_e,g_fname_e,...
                                 theta,phi,U_e,Inf,Inf,options_e,x0_e,fb_e);


%%
%-----------------------------------------------------------------------
%%% SPECIFY PRIORS ON EXTENDED MODEL
%-----------------------------------------------------------------------

% Defining Priors
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim_e.n_phi,1); 
priors.muTheta = zeros(dim_e.n_theta,1);
priors.SigmaPhi = 1e2*eye(dim_e.n_phi);
priors.SigmaTheta = 1e2*eye(dim_e.n_theta);
% Priors on initial 
priors.muX0 = ones(dim_e.n,1)*0;
priors.SigmaX0 = 0e4*eye(dim_e.n);
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;

options_e.priors = priors;
%-----------------------------------------------------------------------
%%% MODEL INVERSION
%-----------------------------------------------------------------------

disp('---------- Invert Data ----------')
[posterior,out] = VBA_NLStateSpaceModel(Y,U, f_fname_e,g_fname_e,dim_e,options_e);
%[posterior,out] = VBA_NLStateSpaceModel(y,u, f_fname_e,g_fname_e,dim_e,options_e);



