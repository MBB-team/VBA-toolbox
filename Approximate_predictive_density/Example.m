% Computes predictive density for a delay discounting model

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



%---- Definition of the model : Hyperbolic delay discount

g_fname = @g_1Dhyp; % observation function, Parameters : [K, beta] => n_phi = 2
dim = struct('n',0,'n_theta',0,'n_phi',2,'p',N,'n_t',1);
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.SigmaPhi = 1e4*eye(dim.n_phi);
% priors.SigmaPhi(end,end) = 0; % Do not infer beta!
% No state noise for deterministic update rules
% Options for inversion
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.dim = dim;
options.verbose = 0;
inG.T = T; % 2 * N (times of reception)
inG.V = OV;  % 2 * N  (objective values)
options.inG = inG; 
K_hyp = 0.1;
beta_hyp = 2;
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;

%-------------- Compute Predictive density

N = 100;
% through MCMC sampling
[pX,gX,pY,gY,X,Y] = get_MCMC_predictiveDensity([],g_fname,u,dim.n_t,options,dim,N);
% through the Laplace approximation
[muy,Vy,iVp] = VBA_getLaplace(u,[],g_fname,dim,options);
