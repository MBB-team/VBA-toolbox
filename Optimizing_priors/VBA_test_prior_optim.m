% Script to test the optimization of priors for model comparison
% In this example, I want to compare two models of learning and decision
% making on an operant learning task (as in Pessiglione 2006)

clear all
close all
clc
 
%---- INDICES IN THE SCRIPT
% - i_m : refers to the index of models
% - i_inv : for each model, inversions to be performed on the parameter
% grid are ordered. i_inv is the index over this list of set of parameters



%---------------------------------------------------------------
% Declaring task contingencies
%---------------------------------------------------------------

Ntrials = 60;
p1_1 = 0.7; % probability of outcome 1 for choice 1
p1_2 = 0.3; % probability of outcome 1 for choice 2

u1 = [rand(1,Ntrials)<p1_1];  % outcome category for alternative 1;
u2 = [rand(1,Ntrials)<p1_2];  % outcome category for alternative 2;

o1 = 1; % value of outcome 1
o2 = 0; % value of outcome 2  

 R = rand(2,Ntrials); % Rewards for each alternatives
 I1 = find(R==1);
 I0 = find(R==0);
 R(I1) = o1;
 R(I0) = o2;
 
%---------------------------------------------------------------
% Declaring models
%---------------------------------------------------------------

Nmodels = 2;
models = cell(1,Nmodels); % containing information about the models
simulations = cell(1,Nmodels); % containing results of simulations per model

%---- MODEL 1 : REINFORCEMENT LEARNING MODEL

%-- Model definition (single session)
f_fname = @f_Qlearn_2Q;
g_fname = @g_softmax_2Q;

dim_output = 1; % the choice made
dim_data = 2; % index of sequence(not used for sim)/ rewards for both alternatives
dim = struct('n',2,...  %( 2 (Qvalues)
    'p',dim_output,... % total output dimension
    'n_theta',1,... % evolution parameters
    'n_phi', 1,... % observation parameters
    'n_t',Ntrials);
options.dim = dim;

options.inF = [];
options.inG = [];
options.binomial = 1; % Dealing with binary data
options.skipf = zeros(1,Ntrials);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

% For any predefined rewards for all alternatives
h_fname = @h_reward_2Q;
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = [2]; % indices where put feedbacks in the experimenter data matrix u

models{1}.options = options;
models{1}.fb = fb;
models{1}.f_fname = f_fname;
models{1}.g_fname = g_fname;

%---- MODEL 2 : PROBABILITY LEARNING MODEL

f_fname = @f_OpLearn_2p; % evolution function for the learning of probability
g_fname = @g_softmax_EU_2p; % softmax decision based on probabilities of outcome 1
% parameters : 1=inverse temperature, 2= asymmetry parameter in prospect
% theory
h_fname = @h_choice_outcome_2p; % feedback is the outcome of the chosen action

dim = struct('n',2*5,... % number of hidden states in the probability learning model
             'p', 1,...
             'n_theta',3,...
             'n_phi',2,...
             'n_t',Ntrials);
         
        
fb.inH.u0 = [u1;u2]; % definition of the binary time-series to be predicted
fb.h_fname = h_fname;
fb.indy = 1; % where to write subject choice in vector u
fb.indfb = 2; % where to write subject choice 

% defining the utility function
u_fname = @u_prospect_theory; % handle of utility function
inG.u_fname = u_fname;
inG.o1 = o1;
inG.o2 = o2;

% simulation parameters % See  Mathys, Daunizeau et al. 2010 for a
% detailled description of parameters
inF.lev2 = 1; % remove 3rd level (volatility learning)
inF.kaub = 1.4;
inF.thub = 1;
inF.rf = -1;
inG.respmod = 'taylor';

options.inF = inF;
options.inG = inG;
options.skipf = zeros(1,Ntrials);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.
       
options.dim = dim;         
options.binomial = 1; % Dealing with binary data
         
models{2}.options = options;
models{2}.fb = fb;
models{2}.f_fname = f_fname;
models{2}.g_fname = g_fname;

%---------------------------------------------------------------
% Running multiple simulations for simulation
%---------------------------------------------------------------

Nsim = [5,5]; % Number of simulations per model

%---- MODEL 1 : REINFORCEMENT LEARNING MODEL

%-- Simulation parameters

alpha = 0.1; % learning rate
beta = 3; % inverse temperature
x0 = [0;0];

theta = sigm(alpha,struct('INV',1)); % transformed parameters
phi = log(beta); % transformed parameters

param.theta = theta*ones(1,Nsim(1)); % might be different parameters for the different simulations
param.phi = phi*ones(1,Nsim(1)); % might be different parameters for the different simulations
param.x0 = x0*ones(1,Nsim(1)); % idem


options = models{1}.options;
fb = models{1}.fb;
f_fname = models{1}.f_fname;
g_fname = models{1}.g_fname;

%-- Running  simulations

Simulations = cell(1,Nsim(1)); % Containing results of simulations

for i_sim = 1 : Nsim(1)
    
    theta = param.theta(:,i_sim);
    phi = param.phi(:,i_sim);
    x0 = param.x0(:,i_sim);
    
    fb.inH.u0 = R; % with reversals
    u = [zeros(2,Ntrials)];
    
    [y,x,x0,eta,e,u] = simulateNLSS_fb(Ntrials,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);    
    isYout = zeros(1,Ntrials);
    
    % saving
    simulations{1}{i_sim}.x0 = x0;
    simulations{1}{i_sim}.theta = theta;
    simulations{1}{i_sim}.phi = phi;
    simulations{1}{i_sim}.u = u;
    simulations{1}{i_sim}.y = y;
    simulations{1}{i_sim}.isYout = isYout;
    
end


%---- MODEL 2 : PROBABILITY LEARNING MODEL
%-- Simulation parameters

 phi = [log(2);... inverse temperature
       1];       % Asymetry in utility of prospect theory
 theta = [1;-4;-1];
 
x0 = repmat([0.5;0;0;1;log(4)],2,1);

param.theta = theta*ones(1,Nsim(1)); % might be different parameters for the different simulations
param.phi = phi*ones(1,Nsim(1)); % might be different parameters for the different simulations
param.x0 = x0*ones(1,Nsim(1)); % idem

options = models{2}.options;
fb = models{2}.fb;
f_fname = models{2}.f_fname;
g_fname = models{2}.g_fname;

%-- Running  simulations

Simulations = cell(1,Nsim(2)); % Containing results of simulations

for i_sim = 1 : Nsim(2)
    
    theta = param.theta(:,i_sim);
    phi = param.phi(:,i_sim);
    x0 = param.x0(:,i_sim);
    
    fb.inH.u0 =  [u1;u2]; % with reversals
    u = [zeros(3,Ntrials)];
    
    [y,x,x0,eta,e,u] = simulateNLSS_fb(Ntrials,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);    
    isYout = zeros(1,Ntrials);
    
    % saving
    simulations{2}{i_sim}.x0 = x0;
    simulations{2}{i_sim}.theta = theta;
    simulations{2}{i_sim}.phi = phi;
    simulations{2}{i_sim}.u = u;
    simulations{2}{i_sim}.y = y;
    simulations{2}{i_sim}.isYout = isYout;
    

end

%---------------------------------------------------------------
% Declaring priors to optimize and a space on which to optimize
%---------------------------------------------------------------

optim_options = struct();

%---------------------------------------------------------------
% Declaring the priors that won't be optimized
%---------------------------------------------------------------

priors = cell(1,Nmodels);

%---- MODEL 1 : REINFORCEMENT LEARNING MODEL
% mean fixed for all priors for the RL model
% std fixed for all priors for the RL model

% By default, all priors are set to dirac.
% Optimized priors will be set to different values during optimization
dim = models{1}.options.dim;
x0 = [0;0];
priors{1}.muX0 = x0;
priors{1}.SigmaX0 = 0e4*eye(dim.n);
priors{1}.muPhi = 0*ones(dim.n_phi,1);
priors{1}.SigmaPhi = 0e4*eye(dim.n_phi);
priors{1}.muTheta = 0*ones(dim.n_theta,1);
priors{1}.SigmaTheta = 0e4*eye(dim.n_theta);

%---- MODEL 2 : PROBABILITY LEARNING MODEL

% By default, all priors are set to dirac on 0
% Optimized priors will be set to different values during optimization
dim = models{2}.options.dim;
x0 = repmat([0.5;0;0;1;log(4)],2,1);
priors{2}.muX0 = x0;
priors{2}.SigmaX0 = 0e4*eye(dim.n);
priors{2}.muPhi = 0*ones(dim.n_phi,1);
priors{2}.SigmaPhi = 0e4*eye(dim.n_phi);
priors{2}.muTheta = 0*ones(dim.n_theta,1);
priors{2}.SigmaTheta = 0e4*eye(dim.n_theta);

optim_options.priors = priors;



%---------------------------------------------------------------
% Declaring the parameters on which the optimization will be performed
%---------------------------------------------------------------

% Here, as an example we want only want to optimize on parameters of the
% two models.


param2optim = cell(1,Nmodels);
% Description of the content of param2optim
% - param2optim{i}.(x/theta/phi).(mu/s).(ind/values/step/bounds)

%---- MODEL 1 : REINFORCEMENT LEARNING MODEL

priors = optim_options.priors{1};
dim = models{1}.options.dim;

prior2optim{1}.theta.mu.ind = [1];
prior2optim{1}.theta.mu.step = 1;
prior2optim{1}.theta.mu.bounds = [-10,10];

prior2optim{1}.theta.s.ind = [1];
prior2optim{1}.theta.s.step = 1;
prior2optim{1}.theta.s.bounds = [-10,10];

prior2optim{1}.phi.mu.ind = [1];
prior2optim{1}.phi.mu.step = 1;
prior2optim{1}.phi.mu.bounds = [-10,10];

prior2optim{1}.phi.s.ind = [1];
prior2optim{1}.phi.s.step = 1;
prior2optim{1}.phi.s.bounds = [-10,10];

%---- MODEL 2 : PROBABILITY LEARNING MODEL

priors = optim_options.priors{2};

prior2optim{2}.theta.mu.ind = [1:3];
prior2optim{2}.theta.mu.step = [1,1,1];
prior2optim{2}.theta.mu.bounds = [-10,10;
                                  -10,10;
                                  -10,10];

prior2optim{2}.theta.s.ind = [1:3];
prior2optim{2}.theta.s.step = [1,1,1];
prior2optim{2}.theta.s.bounds = [-10,10;
                                  -10,10;
                                  -10,10];
                              
prior2optim{2}.phi.mu.ind = [1];
prior2optim{2}.phi.mu.step = 1;
prior2optim{2}.phi.mu.bounds = [-10,10];
                              
prior2optim{2}.phi.s.ind = [1];
prior2optim{2}.phi.s.step = 1;
prior2optim{2}.phi.s.bounds = [-10,10];
 
optim_options.prior2optim = prior2optim;

%---------------------------------------------------------------
% Declaring the optimization criteria
%---------------------------------------------------------------



%---------------------------------------------------------------
% Declaring the optimization method
%---------------------------------------------------------------

optim_options.method = 'grid search';
optim_options.gridsize = 10;

%---------------------------------------------------------------
% Performing the optimization
%---------------------------------------------------------------

results = Optimizing_priors_for_model_comparison(simulations,models,optim_options);
