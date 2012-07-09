% Optimizing priors to compare 2 models of hyperbolic discount
% Cf demo delay
clear all 
close all

%---- Definition of the task contingencies

N = 100; % number of trials
T = zeros(2,N); % time of reception of alternatives
T(1,:) = zeros(1,N); % random choice (uniform)  
T(2,:) = T(1,:) +ones(1,N)*10 + randn(1,N); % T2 = T1 + random delay  
OV = zeros(2,N); % value of alternatives
OV(1,:) = 10; % fixed objective value
OV(2,:) = 20; % fixed objective value
u = [];


% Value now : 10
% value later : 20 discounted

% We want a discount of 1/2 at To = 10
% - hyp : (1+K_hyp*To) = 2 => K_hyp = 1/To = 0.1
% - exp : exp(-K_exp*To) = 2 => K_exp = log(2)/To = 0.07
%---------------- Plotting contingencies

inG.T = T; % 2 * N (times of reception)
inG.V = OV;  % 2 * N  (objective values)
options.inG = inG; 

% Model 1
K_hyp = 0.1;
beta_hyp = 2;
% Model 2
K_exp = 0.07;
beta_exp = 2;



% plotting distribution of the difference of bjective values
p_hyp = g_1Dhyp([],log([K_hyp;beta_hyp]),[],inG)

dov = OV(1,:)-OV(2,:)
dsv = OV(1,:)./(1+K_hyp*T(1,:))-OV(2,:)./(1+K_hyp*T(2,:))

figure
subplot(1,2,1)
hist(dsv)
title('difference of subjective value')
subplot(1,2,2)
hist(p_hyp)
title('proba of choice 1')

% plotting distribution of the difference of bjective values
p_exp = g_1Dexp([],log([K_exp;beta_exp]),[],inG)

dov = OV(1,:)-OV(2,:)
dsv = OV(1,:).*exp(-K_exp.*T(1,:))-OV(2,:).*exp(-K_exp.*T(2,:))

figure
subplot(1,2,1)
hist(dsv)
title('difference of subjective value')
subplot(1,2,2)
hist(p_exp)
title('proba of choice 1')





%%



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
priors.a_sigma = 1;
priors.b_sigma = 1
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
priors.SigmaPhi = 1e4*eye(dim.n_phi);
%priors.SigmaPhi(end,end) = 0; % Do not infer beta!
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;
priors.a_sigma = 1;
priors.b_sigma = 1;

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
[o1,o2,o3] = VBA_check([],u,M{1}.f_fname,M{1}.g_fname,M{1}.options.dim,M{1}.options);
density{1} = o1.priors;
density{1}.a_sigma = 1;
density{1}.b_sigma = 1;

density{1}.muPhi = log([K_hyp;beta_hyp]); 
density{1}.SigmaPhi = 0e1*eye(M{1}.options.dim.n_phi)/10;
density{1}.a_alpha = Inf;
density{1}.b_alpha = 0;

% Model 2
[o1,o2,o3] = VBA_check([],u,M{2}.f_fname,M{2}.g_fname,M{2}.options.dim,M{2}.options);
density{2} = o1.priors;
density{2}.a_sigma = 1;
density{2}.b_sigma = 1;

density{2}.muPhi = log([K_exp;beta_exp]); 
density{2}.SigmaPhi = 0e1*eye(M{2}.options.dim.n_phi);
density{2}.a_alpha = Inf;
density{2}.b_alpha = 0;

% %---------------- Launching the optimization

 partition  = [];
 Nsim = 1;
 n_t = 1; % temporal  dimension in sessions...
 
 m_gen = 1;
 opt = M{m_gen}.options;
 opt.priors = density{m_gen};
 
 
[pX,gX,pY,gY,X,Y] = get_MCMC_predictiveDensity(M{m_gen}.f_fname,...
                                               M{m_gen}.g_fname,...
                                               u,...
                                               n_t,...
                                               opt,...
                                               M{m_gen}.options.dim,...
                                               Nsim);

                                           i_sim= 1
                                           m_inv = 1;
                                           
       opt = M{m_gen}.options;
                                     
 [posterior,out] = VBA_NLStateSpaceModel(Y(:,:,i_sim),...
                                                        u,...
                                                        M{m_inv}.f_fname,...
                                                        M{m_inv}.g_fname,...
                                                        M{m_inv}.options.dim,...
                                                        opt);
