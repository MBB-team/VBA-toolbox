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


%--------------------------------------------------
%
% EV










clear all 
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATA SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulate data from single session
Ntrials = 10;
% Simulation parameters
alpha = 0.2; % learning rate
beta =10; % inverse temperature
Q0 = [0;0];

% Probabilist ic reward. Correlated probabilities (complementary p1 = 1-p2)
% Fixed probabilities (p1 =0.8) with probability reversal at mid course
R = Inf*zeros(2,Ntrials); % Rewards for each alternatives
R(1,1:Ntrials) = 1;
R(2,1:Ntrials) = 0;


% simpler R
%R(1,:) = 1;
%R(2,:) = 0;

x0 = [0;0];
n_t = Ntrials;
% allocate feedback struture for simulations

% For any predefined rewards for all alternatives
h_fname = @h_reward_2Q;
fb.inH.u0 = R;

fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = [2]; % indices where put feedbacks in the experimenter data matrix u

f_fname = @f_Qlearn_2Q;
g_fname = @g_softmax_2Q;
theta = sigm(alpha,struct('INV',1));
phi = [log(beta);0.5];
u = [nan*zeros(2,Ntrials)];

alpha = Inf;
sigma =  Inf;

dim_output = 1; % the choice
dim = struct('n',2,...  %( 2 (Qvalues) 
             'p',dim_output,... % total output dimension
             'n_theta',1,... % evolution parameters
             'n_phi', 2,... % observation parameters
             'n_t',Ntrials);


options.inF.param_transform.type = 'modified sigmoid';
options.inF.param_transform.a = 0;
options.inF.param_transform.b = 1;

options.inG.param_transform.type = 'modified sigmoid';
options.inG.param_transform.a = 1;
options.inG.param_transform.b = 5;

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
priors.SigmaPhi = 1e0*eye(dim.n_phi);
priors.SigmaTheta = 1e0*eye(dim.n_theta);
% priors.SigmaPhi(2,2) = 0;
% Priors on initial 
priors.muX0 = ones(dim.n,1)*0;
priors.SigmaX0 = 0e0*eye(dim.n);
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;

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
