%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this script, I inverse a simple reinforcement learning model with
% softmax decision rule.
% The task consists of a sequential decision task among two alternatives
% that are probabilistically rewarded 
% p(reward) = 0.8 for the first // p(reward) = 0.2 for the second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all
clc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATA SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate data from multiple sessions
Nsessions =4;

%--------------------------------- Shared by all sessions for simulation

Ntrials = 60;
% Simulation parameters
alpha = 0.1; % learning rate
beta =3; % inverse temperature
U = []; Y = []; IsYout = [];
x0 = [0;0];

f_fname = @f_Qlearn_2Q;
g_fname = @g_softmax_2Q;
theta = sigm(alpha,struct('INV',1));
phi = log(beta);
alpha = Inf;
sigma =  Inf;

dim_output = 1; % the delta reaction-time
dim_data = 2; % index of sequence(not used for sim)/ rewards for both alternatives
dim = struct('n',2,...  %( 2 (Qvalues) * Nsessions
             'p',dim_output,... % total output dimension
             'n_theta',1,... % evolution parameters
             'n_phi', 1,... % observation parameters
             'n_t',Ntrials);
options.dim = dim;

options.inF = [];
options.inG = [];
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.isYout = zeros(1,Ntrials); % Excluding data points
options.dim = dim;
                
options.skipf = zeros(1,Ntrials);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.


% For any predefined rewards for all alternatives
h_fname = @h_reward_2Q;
u0 = [ones(1,Ntrials/3)]; % possible feedbacks
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = [2]; % indices where put feedbacks in the experimenter data matrix u


%-----------------------------------------------------------------------

for i = 1 : Nsessions

% Probabilistic reward. Correlated probabilities (complementary p1 = 1-p2)
% Fixed probabilities (p1 =0.8) with probability reversal at mid course
R = rand(2,Ntrials); % Rewards for each alternatives
R(1,1:Ntrials/2) = R(1,1:Ntrials/2)<0.8;R(1,Ntrials/2+1:end) = R(1,Ntrials/2+1:end)<0.2;
R(2,1:Ntrials/2) = R(2,1:Ntrials/2)<0.2;R(2,Ntrials/2+1:end) = R(2,Ntrials/2+1:end)<0.8;

fb.inH.u0 = R; % with reversals
u = [zeros(2,Ntrials)];

[y,x,x0,eta,e,u] = simulateNLSS_fb(Ntrials,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);

isYout = zeros(1,Ntrials);
U = [U;u];
Y = [Y;y];
IsYout = [IsYout;isYout];


end


%%
%-----------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPECIFY MODEL FOR MULTIPLE SESSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dim_output = 1; % the delta reaction-time
dim_data = 2; % reward/index of sequence 
dim = struct('n',2*Nsessions,...  %( 2 (special for RL) * 2 (Qvalues) * Nsessions
             'n_ps', 2,... % number of hidden state per session
             'p',Nsessions*dim_output,... % total output dimension
             'p_ps',dim_output,... % output dimension per sessions
             'n_theta',1,... % evolution parameters
             'n_phi', 1,... % observation parameters
             'n_data_ps',dim_data,... % data dimension per session
             'n_t',Ntrials,...
             'n_sess', Nsessions); %

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
in_sessions = struct();
in_sessions.f_fname = @f_Qlearn_2Q;
in_sessions.g_fname = @g_softmax_2Q;
in_sessions.dim = dim;
in_sessions.ind.theta =  ones(Nsessions,1)*[1];
in_sessions.ind.phi =   ones(Nsessions,1)*[1];
in_sessions.binomial = 1;
in_sessions.inF = [];
in_sessions.inG = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL INVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[posterior,out] = inversion_multiple_sessions(in_sessions, Y, U, IsYout, priors);


