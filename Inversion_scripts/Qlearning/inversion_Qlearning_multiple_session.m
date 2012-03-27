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
Ntrials = 200;
% Simulation parameters
alpha = 0.1; % learning rate
beta =3; % inverse temperature
Q0 = [0;0];
Q = [];
U = []; Y = []; IsYout = [];
for i = 1 : Nsessions
    
R = rand(2,Ntrials); % Rewards for each alternatives
R(1,1:Ntrials/2) = R(1,1:Ntrials/2)<0.8;R(1,Ntrials/2+1:end) = R(1,Ntrials/2+1:end)<0.2;
R(2,1:Ntrials/2) = R(2,1:Ntrials/2)<0.2;R(2,Ntrials/2+1:end) = R(2,Ntrials/2+1:end)<0.8;

[A,r,q] = simulateQLearning(Q0,alpha,beta,R); % simulate RL model
u = [A;r];
y = A;
isYout = zeros(1,Ntrials);
U = [U;u];
Y = [Y;y];
IsYout = [IsYout;isYout];
Q = [Q;q;q];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPECIFY MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dim_output = 1; % the delta reaction-time
dim_data = 2; % reward/index of sequence 
dim = struct('n',2*2*Nsessions,...  %( 2 (special for RL) * 2 (Qvalues) * Nsessions
             'n_ps', 2*2,... % number of hidden state per session
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL INVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[posterior,out] = inversion_multiple_sessions(in_sessions, Y, U, IsYout, priors);
sigm(posterior.muTheta,struct('INV',0))
log(posterior.muPhi)




