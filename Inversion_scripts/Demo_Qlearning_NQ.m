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
Ntrials = 300;
Na = 3; % Case of 3 Qvalues

% Simulation parameters
alpha = 0.1; % learning rate
beta =5; % inverse temperature
Q0 = [0;0;0]; % Qvalues initialy set to zeros

% Probabilistic reward.
% 3 blocks. One option alternatively highly rewarded ([0.8,0.2,0.2] 
I = ones(1,Ntrials/3);
R = rand(3,Ntrials); % Rewards for each alternatives
R = R<[0.8*I,0.2*I,0.2*I;
       0.2*I,0.8*I,0.2*I;
       0.2*I,0.2*I,0.8*I];


x0 = [0;0;0;0;0;0];
n_t = Ntrials;


% For any predefined rewards for all alternatives
h_fname = @h_reward_NQ;
u0 = [ones(1,Ntrials/3)]; % possible feedbacks
fb.inH.u0 = R; % with reversals
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = [2,3,4]; % indices where put feedbacks in the experimenter data matrix u
fb.Na = Na;

f_fname = @f_Qlearn_NQ;
g_fname = @g_softmax_NQ;
theta = sigm(alpha,struct('INV',1));
phi = log(beta);
u = [zeros(1,Ntrials);R];
u(1,1)=1; % This forces the first decision
alpha = Inf;
sigma =  Inf;

dim_output = 1; % the delta reaction-time
dim_data = 4; % index of sequence(not used for sim)/ rewards for both alternatives
dim = struct('n',Na*2,...  %( 2 (special for RL) * 2 (Qvalues) * Nsessions
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
options.ind_copy = [1,3,5;...
                    2,4,6]; % copy 1 to 2, copy 3 to 4

                
 [y,x,x0,eta,e,u] = simulateNLSS_fb_multinomial(n_t,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);
%[y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);
%[y,x,x0,eta,e,u] = simulate_model_multinomial(n_t,Na,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

figure
hold on
plot(x(1,:),'b')
plot(x(3,:),'r')
plot(x(5,:),'g')
plot(y,'kx')
legend('0','1')
%%


% Defining Priors
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = phi;%zeros(dim.n_phi,1); 
priors.muTheta = theta;%zeros(dim.n_theta,1);
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

g_fname = @g_softmax_NQ_inv;


%%%% Options for inversion
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.priors = priors;

y = ones(1,Ntrials);
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)

