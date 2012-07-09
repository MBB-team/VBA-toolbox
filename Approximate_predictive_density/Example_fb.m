% Computes predictive density of a reinforcement learning model
% See Demo_QLearning_2Q.m

clear all 
close all


%---- Definition of the task contingencies
% Simulate data from single session
Ntrials = 100;
% Simulation parameters
alpha = 0.2; % learning rate
beta =1; % inverse temperature
Q0 = [0;0];

% Probabilist ic reward. Correlated probabilities (complementary p1 = 1-p2)
% Fixed probabilities (p1 =0.8) with probability reversal at mid course
R = rand(2,Ntrials); % Rewards for each alternatives
R(1,1:Ntrials/2) = R(1,1:Ntrials/2)<0.8;R(1,Ntrials/2+1:end) = R(1,Ntrials/2+1:end)<0.2;
R(2,1:Ntrials/2) = R(2,1:Ntrials/2)<0.2;R(2,Ntrials/2+1:end) = R(2,Ntrials/2+1:end)<0.8;

x0 = [0;0];
n_t = Ntrials;
% allocate feedback struture for simulations

% For any predefined rewards for all alternatives
h_fname = @h_reward_2Q;
u0 = [ones(1,Ntrials/3)]; % possible feedbacks
fb.inH.u0 = R;
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = [2]; % indices where put feedbacks in the experimenter data matrix u

%---- Definition of the model : Reinforcement learning


f_fname = @f_Qlearn_2Q;
g_fname = @g_softmax_2Q;

u = [zeros(2,Ntrials)];

dim_output = 1; % the delta reaction-time
dim_data = 3; % index of sequence(not used for sim)/ rewards for both alternatives
dim = struct('n',2,...  %( 2 (Qvalues) * Nsessions
             'p',dim_output,... % total output dimension
             'n_theta',1,... % evolution parameters
             'n_phi', 1,... % observation parameters
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

% Defining parameters for sessions
options.binomial = 1; % Dealing with binary data
options.priors = priors;
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;

%-------------- Compute Predictive density
%%
% through MCMC sampling
N = 100; % number of samples
[pX,gX,pY,gY,X,Y,U] = get_MCMC_predictiveDensity_fb(f_fname,g_fname,u,dim.n_t,options,dim,fb,N);
