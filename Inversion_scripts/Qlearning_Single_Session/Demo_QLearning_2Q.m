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
Ntrials = 500;
% Simulation parameters
alpha = 0.05; % learning rate
beta =3; % inverse temperature
Q0 = [0;0];

% Probabilistic reward. Correlated probabilities (complementary p1 = 1-p2)
% Fixed probabilities (p1 =0.8) with probability reversal at mid course
R = rand(2,Ntrials); % Rewards for each alternatives
R(1,1:Ntrials/2) = R(1,1:Ntrials/2)<0.8;R(1,Ntrials/2+1:end) = R(1,Ntrials/2+1:end)<0.2;
R(2,1:Ntrials/2) = R(2,1:Ntrials/2)<0.2;R(2,Ntrials/2+1:end) = R(2,Ntrials/2+1:end)<0.8;

[A,r,q] = simulate_QLearning_2Q(Q0,alpha,beta,R); % simulate RL model
u = [A;r];
y = A;
isYout = zeros(1,Ntrials);

plot(q')
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPECIFY MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim_output = 1; % the delta reaction-time
dim_data = 2; % reward/index of sequence 
dim = struct('n',2*2,...  %( 2 (special for RL) * 2 (Qvalues) * Nsessions
             'p',dim_output,... % total output dimension
             'n_theta',1,... % evolution parameters
             'n_phi', 1,... % observation parameters
             'n_t',Ntrials);

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
options.inF = [];
options.inG = [];
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.isYout = isYout; % Excluding data points


theta = sigm(alpha,struct('INV',1));
phi = log(beta);
[posterior,out] = VBA_NLStateSpaceModel(y,u,@f_Qlearn_2Q,@g_softmax_2Q,dim,options);

% Reformating for plot (getting rid of copy of Qvalues)
posterior2 = posterior;
posterior2.muX = posterior.muX([1,3],:);
for i = 1:Ntrials
posterior2.SigmaX.current{i} = posterior.SigmaX.current{i}([1,3],[1,3]);
end
displayResults(posterior2,out,y,[q(1,:);q(2,:)],[0;0],theta,phi,Inf,Inf)


