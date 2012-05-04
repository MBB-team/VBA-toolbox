%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this script, inversion of a reinforcement learning model with
% softmax decision rule.
% RL model : 2 Qvalues, feedback for the selected action, Rescorla Wagner
% update rule
% Task : sequential decision task among two alternatives
% that are probabilistically rewarded 
% p(reward) = 0.8 for the first // p(reward) = 0.2 for the second
% Probability reversal at mid course
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATA SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulate data from single session
Ntrials = 90;
% Simulation parameters
alpha = 0.2; % learning rate
beta =3; % inverse temperature
Q0 = [0;0];

% Probabilistic reward. Correlated probabilities (complementary p1 = 1-p2)
% Fixed probabilities (p1 =0.8) with probability reversal at mid course
R = rand(2,Ntrials); % Rewards for each alternatives
R(1,1:Ntrials/2) = R(1,1:Ntrials/2)<0.8;R(1,Ntrials/2+1:end) = R(1,Ntrials/2+1:end)<0.2;
R(2,1:Ntrials/2) = R(2,1:Ntrials/2)<0.2;R(2,Ntrials/2+1:end) = R(2,Ntrials/2+1:end)<0.8;


% simpler R
%R(1,:) = 1;
%R(2,:) = 0;

x0 = [0;0];
n_t = Ntrials;

% allocate feedback struture for simulations

%{ 
%This is for coupled rewards!
h_fname = @h_truefalse;
u0 = [ones(1,50)]; % possible feedbacks
fb.inH.u0 = [u0,~u0,u0,~u0,u0,~u0]; % with reversals
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = 2;
%}

% For any predefined rewards for all alternatives
h_fname = @h_reward_2Q;
u0 = [ones(1,Ntrials/3)]; % possible feedbacks
fb.inH.u0 = [u0,~u0,u0;
            ~u0,u0,~u0]; % with reversals
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = [2,3]; % indices where put feedbacks in the experimenter data matrix u

%%% Comparison to what I had done
% feedbacks where specified in advance in u. And simulation just added the
% the generated data in the first line.
% here feedbacks are specified in the feedback structure




f_fname = @f_Qlearn_2Q;
g_fname = @g_softmax_2Q;
theta = sigm(alpha,struct('INV',1));
phi = log(beta);
u = [zeros(1,Ntrials);R];
alpha = Inf;
sigma =  Inf;

dim_output = 1; % the delta reaction-time
dim_data = 3; % index of sequence(not used for sim)/ rewards for both alternatives
dim = struct('n',2,...  %( 2 (Qvalues) * Nsessions
             'p',dim_output,... % total output dimension
             'n_theta',1,... % evolution parameters
             'n_phi', 1,... % observation parameters
             'n_t',Ntrials);
         
options.inF = [];
options.inG = [];
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.isYout = zeros(1,Ntrials); % Excluding data points
options.dim = dim;
                
options.skipf = zeros(1,n_t);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

[y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);


plot(x')
%%

% Defining Priors
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaPhi = 1e2*eye(dim.n_phi);
priors.SigmaTheta = 1e2*eye(dim.n_theta);
% Priors on initial 
priors.muX0 = ones(dim.n,1)*0;
priors.SigmaX0 = 0e4*eye(dim.n);
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;
% Defining parameters for sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL INVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Options for inversion
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.priors = priors;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)
