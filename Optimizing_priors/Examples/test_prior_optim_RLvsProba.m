% Optimizing priors to compare 2 models of operant learning

clear all 
close all

%-------------------------------------------------
%---- Definition of the task contingencies
%-------------------------------------------------

% repeated choice among two alternatives, each leading to the same two
% possible outcomes but with different probabilities
u1 = [rand(1,50)<0.8,rand(1,50)<0.2];  % outcome category for alternative 1; (for proba tracking)
u2 = [rand(1,50)<0.2,rand(1,50)<0.8];  % outcome category for alternative 2;
u = [u1;u2];
Ntrials = length(u1);
% Definition of the outcome values for each of the two possible outcomes
o1 = 1;
o2 = -1;
R = zeros(size(u));
R(find(u==0)) = o1; % outcome value for alternative 1 (for RL)
R(find(u==1)) = o2; % outcome value for alternative 2

%-------------------------------------------------
%---- Definition of the models
%-------------------------------------------------

M = cell(1,2);

%---- Model 1 : Reinforcement learning

% allocate feedback struture for simulations
% For any predefined rewards for all alternatives
h_fname = @h_reward_2Q;
fb.inH.u0 = R;
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = [2]; % indices where put feedbacks in the experimenter data matrix u

M{1}.g_fname = @g_softmax_2Q;
M{1}.f_fname = @f_Qlearn_2Q;
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
                
options.skipf = zeros(1,Ntrials);
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
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;
M{1}.options = options;

%---- Model 2 : Probability learning

% evolution, observation and feedback functions
M{2}.f_fname = @f_OpLearn_2p; % evolution function for the learning of probability
M{2}.g_fname = @g_softmax_EU_2p; % softmax decision based on probabilities of outcome 1
% parameters of g: 1=inverse temperature, 2= asymmetry parameter in prospect
% theory

% defining the utility function
inG.u_fname = @u_prospect_theory;% handle of utility function
inG.o1 = o1;
inG.o2 = o2;

% simulation parameters
inF.lev2 = 0; % remove 3rd level (volatility learning)
inF.kaub = 1.4; % upper bound for Kappa
inF.thub = 1; % upper bound for Theta
inF.rf = -1;
inG.respmod = 'taylor';

% choose initial conditions
x0 = repmat([0.5;0;0;1;log(4)],2,1); % identical initial conditions for both models

u = [zeros(2,Ntrials)]; % vector u consists in choice, and outcome category
dim = struct('n',2*5,... % number of hidden states (5 by choice)
    'p',1,... % total output dimension
    'n_theta',3,... % number of evolution parameters
    'n_phi',2,... % number of observation parameters
    'n_t',Ntrials); % number of time trials


options.binomial = 1; % we here deal with binary choices
options.inF = inF;
options.inG = inG;
options.dim = dim;
options.skipf = zeros(1,length(u));
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

% Defining Priors
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1);
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaPhi = 1e2*eye(dim.n_phi);
priors.SigmaTheta = 1e2*eye(dim.n_theta);
% Priors on initial
priors.muX0 = x0;%ones(dim.n,1)*0;
priors.SigmaX0 = 0e4*eye(dim.n);
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;
%%%% Options for inversion
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data

options.priors = priors;
M{2}.options = options;


%%
%----------------- Density for data simulation

density = cell(1,2);

%---- Model 1 : Reinforcement learning
% For the simulations of a RL model, we choose a rather fast learner but
% also allow for slower ones by setting a large variance
% to have an idea of the prior on alpha : X = sigm(-0.8473 + randn(1,1e6));hist(X,100)
% to have an idea of the prior on beta : X = exp(0 + randn(1,1e6));hist(X,1000)

alpha = 0.3; % learning rate
beta = 1; % noise level

density{1} = M{1}.options.priors; % loading standard priors
density{1}.muPhi = [log(beta)];
density{1}.SigmaPhi = 1e1*eye(M{1}.options.dim.n_phi);
density{1}.muTheta = [sigm(alpha,struct('INV',1))];
density{1}.SigmaTheta = 1e1*eye(M{1}.options.dim.n_theta);

density{1}.muX0 = ones(dim.n,1)*(o1+o2)/2; % prior Qvalues :mean of outcomes
density{1}.SigmaX0 = ones(dim.n,1)*(o1+o2)/2; % prior Qvalues :mean of outcomes
density{1}.a_sigma = 1;
density{1}.b_sigma = 1;
density{1}.a_alpha = Inf;
density{1}.b_alpha = 0;

%---- Model 2 : Probability learning
% Here, you need to have a good understanding of the model to propose
% parameters leading to an interesting behavior.

% Watch out for the mapping of the parameters
theta = [isgm(1,M{2}.options.inF.kaub);...
        -4;... no mapping for omega
        isgm(-1,M{2}.options.inF.kaub)]; % theta : kappa,omega,theta
phi = [log(3);log(1)]; % inverse temperature in softmax, parameter of utility function

density{2} = M{2}.options.priors; % loading standard priors
density{2}.a_sigma = 1;
density{2}.b_sigma = 1;
density{2}.muPhi = phi;
density{2}.SigmaPhi = 1e1*eye(M{1}.options.dim.n_phi);
density{2}.muTheta = theta;
density{2}.SigmaTheta = 1e1*eye(M{1}.options.dim.n_theta);
density{2}.a_alpha = Inf;
density{2}.b_alpha = 0;

%---------------- Declaring the parameters on which to optimize

% Model 1 : optimize on both observation parameters

priors2optim{1}.phi.mu.ind = [1,2];
priors2optim{1}.phi.mu.step = [1,1];
priors2optim{1}.phi.mu.bounds = [-3,0;0.2,2];

priors2optim{1}.phi.s.ind = [];
priors2optim{1}.phi.s.step = [2,2];
priors2optim{1}.phi.s.bounds = [10,100;10,100];

% Model 2 : optimize on both observation parameters

priors2optim{2}.phi.mu.ind = [1,2];
priors2optim{2}.phi.mu.step = [2,2];
priors2optim{2}.phi.mu.bounds = [-3,0;0.2,2];

priors2optim{2}.phi.s.ind = [];
priors2optim{2}.phi.s.step = [2,2];
priors2optim{2}.phi.s.bounds = [10,100;10,100];

%---------------- Launching the optimization
 partition  = [];
 Nsim = 10;
 
 %%
[optim_priors] = VBA_optimPriors(M,u,partition,density, priors2optim ,Nsim);



