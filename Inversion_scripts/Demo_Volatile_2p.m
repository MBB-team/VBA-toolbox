%----------------------------------------------------------------------
% In this script, inversion of a model learning the binary probability of two cues from feedback
% and selecting action in order to soft-maximize the expected outcome
% utility
% Task : sequential decision task among two alternatives
% that lead to binary probabilistic outcomes (outcome 1 or outcome 2)
%
% -- Remark : 
% The model knows that there are two possible outcomes and knows their
% values
% It learns the probability of each of these two outcomes from the number
% of occurences of the two different outcomes.
% During learning, only the category/index of the outcome matters.
% Only at the time of decision are the actual values of outcomes taken into
% account
%
% -- Reference :
% The model is based on the model used in :
% "A Bayesian foundation for individual learning under uncertainty"
% Christoph Mathys, Jean Daunizeau
%  doi:  10.3389/fnhum.2011.00039
%----------------------------------------------------------------------
clear all
close all
clc

%----------------------------------------------------------------------
% DATA SIMULATION (for a single session)
%----------------------------------------------------------------------

% ------------- Task description & contingencies
% repeated choice among two alternatives, each leading to the same two
% possible outcomes but with different probabilities
u1 = [rand(1,50)<0.8,rand(1,50)<0.2];  % outcome category for alternative 1;
u2 = [rand(1,50)<0.2,rand(1,50)<0.8];  % outcome category for alternative 2;
Ntrials = length(u1);

% allocate feedback struture for simulations
fb.inH.u0 = [u1;u2]; % definition of the binary time-series to be predicted
fb.h_fname = @h_choice_outcome_cat_2p; % feedback is the outcome category of the chosen action
fb.indy = 1; % where to write binary subject choice in vector u
fb.indfb = 2; % where to write binary outcome category in vector u

% Definition of the outcome values for each of the two possible outcomes
o1 = 1;
o2 = -1;

% ------------- Description of the model structure
% evolution, observation and feedback functions
f_fname = @f_OpLearn_2p; % evolution function for the learning of probability
g_fname = @g_softmax_EU_2p; % softmax decision based on probabilities of outcome 1
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
% ---- description of hidden states (for more details see reference)
% - x0(1) : the choice category 
% - x0(2) : mean of belief on probability of output o1 (through sigmoid mapping)
% - x0(3) : variance of belief on probability of output o1 (through sigmoid mapping)
% - x0(4) : mean of belief on volatility
% - x0(5) : variance of belief on volatily
% ---- Explanation of the choices here
% - x0(1) : this has no influence here, first decision is based on x0(2)
% - x0(2) : mapped through a sigmoid, this should reflect the probability of the action to lead to outcome 1.
% Set to 0, means p(o1|a)=p(o2|a)=0.5. Different values could capture an
% initial bias
% - x0(3) : This value (variance) captures the amount of confidence in the prior belief
% on outcome probabilities. A high variance means little confidence. The
% lower the variance, the more the model will rely on its prior on outcome
% probabilities.
% - x0(4) : mean of belief on volatility : the higher the value, the higher
% the believed value of volatility in the environment
% - x0(5) : variance of belief on volatily : amount of confidence on the
% former belief

u = [zeros(2,Ntrials)]; % vector u consists in choice, and outcome category
dim = struct('n',2*5,... % number of hidden states (5 by choice)
    'p',1,... % total output dimension
    'n_theta',3,... % number of evolution parameters
    'n_phi',2,... % number of observation parameters
    'n_t',Ntrials); % number of time trials

theta = [1;-4;-1]; % theta : kappa,omega,theta
phi = [log(3);log(1)]; % inverse temperature in softmax, parameter of utility function

options.binomial = 1; % we here deal with binary choices
options.inF = inF;
options.inG = inG;
options.dim = dim;
options.skipf = zeros(1,length(u));
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

[y,x,x0,eta,e,u] = simulateNLSS_fb(length(u),f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);


% --------- plotting simulation results
figure(1)
hold on
plot(y,'*') % actual binary choices
plot(sigm(x(2,:)),'r') % probability of outcome 1 for choice 1 
plot(sigm(x(7,:)),'b') % probability of outcome 1 for choice 2 
legend('choices','prob outcome1 for cue1','prob outcome1 for cue2')
title('Simulation result')
xlabel('time')


%%
%----------------------------------------------------------------------
% MODEL INVERSION
%----------------------------------------------------------------------

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

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)
