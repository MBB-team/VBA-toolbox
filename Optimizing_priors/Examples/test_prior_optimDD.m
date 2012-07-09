% Optimizing priors to compare 2 models of hyperbolic discount
% Cf demo delay
clear all 
close all
%---- Definition of the task contingencies
% We want a discount of 1/2 at To = 10
% - hyp : (1+K_hyp*To) = 2 => K_hyp = 1/To = 0.1
% - exp : exp(-K_exp*To) = 2 => K_exp = log(2)/To = 0.07

N = 100; % number of trials
T = zeros(2,N); % time of reception of alternatives
T(1,:) = zeros(1,N); % random choice (uniform)  
T(2,:) = T(1,:) +10*ones(1,N) +randn(1,N); % T2 = T1 + random delay  
OV = zeros(2,N); % value of alternatives
OV(1,:) = 10; % fixed objective value
OV(2,:) = 20; % fixed objective value
u = [];
%---- Definition of the models
M = cell(1,2);

%---- Model 1 % Hyperbolic

g_fname = @g_1Dhyp; % observation function, Parameters : [K, beta] => n_phi = 2
dim = struct('n',0,'n_theta',0,'n_phi',2,'p',N,'n_t',1);
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.SigmaPhi = 1e4*eye(dim.n_phi);
%priors.SigmaPhi(end,end) = 0; % Do not infer beta!
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;
% Options for inversion
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.dim = dim;
options.verbose = 0;
inG.T = T; % 2 * N (times of reception)
inG.V = OV;  % 2 * N  (objective values)
options.inG = inG; 

M{1}.options = options;
M{1}.g_fname = g_fname;
M{1}.f_fname = [];

% Model 2

g_fname = @g_1Dexp; % observation function, Parameters : [K, beta] => n_phi = 2
dim = struct('n',0,'n_theta',0,'n_phi',2,'p',N,'n_t',1);
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.SigmaPhi = 1e2*eye(dim.n_phi);
%priors.SigmaPhi(end,end) = 0; % Do not infer beta!
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;
% Options for inversion
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.dim = dim;
options.verbose = 0;
inG.T = T; % 2 * N (times of reception)
inG.V = OV;  % 2 * N  (objective values)
options.inG = inG; 

M{2}.options = options;
M{2}.g_fname = g_fname;
M{2}.f_fname = [];


%----------------- Density for data simulation
density = cell(1,2);

% Model 1
K_hyp = 0.1;
beta_hyp = 2;
[o1,o2,o3] = VBA_check([],u,M{1}.f_fname,M{1}.g_fname,M{1}.options.dim,M{1}.options);
density{1} = o1.priors;
density{1}.a_sigma = 1;
density{1}.b_sigma = 1;
density{1}.muPhi = log([K_hyp;beta_hyp]); 
density{1}.SigmaPhi = 0e1*eye(M{1}.options.dim.n_phi);
density{1}.a_alpha = Inf;
density{1}.b_alpha = 0;
% Model 2
K_exp = 0.07;
beta_exp = 2;
[o1,o2,o3] = VBA_check([],u,M{2}.f_fname,M{2}.g_fname,M{2}.options.dim,M{2}.options);
density{2} = o1.priors;
density{2}.a_sigma = 1;
density{2}.b_sigma = 1;
density{2}.muPhi = log([K_exp;beta_exp]); 
density{2}.SigmaPhi = 0e1*eye(M{2}.options.dim.n_phi);
density{2}.a_alpha = Inf;
density{2}.b_alpha = 0;

%---------------- Declaring the parameters on which to optimize

% Model 1 : optimize on both observation parameters

priors2optim{1}.phi.mu.ind = [1,2];
priors2optim{1}.phi.mu.step = [1,1];
priors2optim{1}.phi.mu.bounds = [-3,0;0.2,2];
% 
% priors2optim{1}.phi.s.ind = [];
% priors2optim{1}.phi.s.step = [2,2];
% priors2optim{1}.phi.s.bounds = [10,100;10,100];
% 
% Model 2 : optimize on both observation parameters

% priors2optim{2}.phi.mu.ind = [1,2];
% priors2optim{2}.phi.mu.step = [2,2];
% priors2optim{2}.phi.mu.bounds = [-3,0;0.2,2];

% priors2optim{2}.phi.s.ind = [];
% priors2optim{2}.phi.s.step = [2,2];
% priors2optim{2}.phi.s.bounds = [10,100;10,100];

%---------------- Launching the optimization
 partition  = [];
 Nsim = 10;
 

 
 %%
[optim_priors] = VBA_optimPriors(M,u,partition,density, priors2optim ,Nsim);



