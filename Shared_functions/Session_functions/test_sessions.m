%-----------------------------------------------------------------------
% In this script, I inverse a simple reinforcement learning model with
% softmax decision rule.
% The task consists of a sequential decision task among two alternatives
% that are probabilistically rewarded 
% p(reward) = 0.8 for the first // p(reward) = 0.2 for the second
%-----------------------------------------------------------------------

clear all 
close all
clc



%-----------------------------------------------------------------------
%%% DATA SIMULATION
%-----------------------------------------------------------------------

% Simulate data from multiple sessions
Nsessions =4;
%--------------------------------- Shared by all sessions for simulation

Ntrials = 60;
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

dim = struct('n',2,...  %( 2 (Qvalues) * Nsessions
             'p',1,... % total output dimension
             'n_theta',1,... % evolution parameters
             'n_phi', 1,... % observation parameters
             'n_t',Ntrials,... % number of trials
             'u',2); % input data
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

%-----------------------------------------------------------------------
%%% SPECIFY MODEL FOR MULTIPLE SESSIONS
%-----------------------------------------------------------------------

dim_e = struct('n_theta',1,... % evolution parameters
             'n_phi', 1);%       observation parameters
             

% Defining parameters for sessions
in_sessions = struct();
in_sessions.n_sess = Nsessions;
in_sessions.f_fname = @f_Qlearn_2Q;
in_sessions.g_fname = @g_softmax_2Q;
in_sessions.dim_e = dim_e;
in_sessions.ind.theta =  ones(Nsessions,1)*[1];
in_sessions.ind.phi =   ones(Nsessions,1)*[1];
in_sessions.binomial = 1;
in_sessions.inF = [];
in_sessions.inG = [];

options.fb = fb;

[ f_fname_e,g_fname_e,dim_e,options_e ] = getModelforSessions( f_fname,g_fname,dim,options,in_sessions);

options_e.isYout = IsYout;
options_e.dim = dim_e;

%%
%-----------------------------------------------------------------------
%%% SIMULATE EXTENDED MODEL
%-----------------------------------------------------------------------
x0_e = zeros(dim_e.n,1);
fb_e = options_e.fb;
U_e = zeros(dim_e.u,Ntrials);

R(1,1:Ntrials/2) = R(1,1:Ntrials/2)<0.8;R(1,Ntrials/2+1:end) = R(1,Ntrials/2+1:end)<0.2;
R(2,1:Ntrials/2) = R(2,1:Ntrials/2)<0.2;R(2,Ntrials/2+1:end) = R(2,Ntrials/2+1:end)<0.8;

fb.inH.u0 = repmat(R,Nsessions,1); % with reversals


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


[posterior,out] = VBA_NLStateSpaceModel(Y,U, f_fname_e,g_fname_e,dim_e,options_e);



